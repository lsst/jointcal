# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest
from unittest import mock

import numpy as np

import lsst.log
import lsst.utils

import lsst.afw.table
import lsst.daf.persistence
from lsst.daf.base import DateTime
import lsst.geom
from lsst.meas.algorithms import getRefFluxField, LoadIndexedReferenceObjectsTask, DatasetConfig
import lsst.pipe.base
import lsst.jointcal
from lsst.jointcal import MinimizeResult
import lsst.jointcal.chi2
import lsst.jointcal.testUtils


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


def make_fake_refcat(center, flux, filterName):
    """Make a fake reference catalog."""
    schema = LoadIndexedReferenceObjectsTask.makeMinimalSchema([filterName])
    catalog = lsst.afw.table.SimpleCatalog(schema)
    record = catalog.addNew()
    record.setCoord(center)
    record[filterName + '_flux'] = flux
    record[filterName + '_fluxErr'] = flux*0.1
    return catalog


class JointcalTestBase:
    def setUp(self):
        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100)
        self.ccdImageList = struct.ccdImageList
        # so that countStars() returns nonzero results
        for ccdImage in self.ccdImageList:
            ccdImage.resetCatalogForFit()

        self.goodChi2 = lsst.jointcal.chi2.Chi2Statistic()
        # chi2/ndof == 2.0 should be non-bad
        self.goodChi2.chi2 = 200.0
        self.goodChi2.ndof = 100

        self.badChi2 = lsst.jointcal.chi2.Chi2Statistic()
        self.badChi2.chi2 = 600.0
        self.badChi2.ndof = 100

        self.nanChi2 = lsst.jointcal.chi2.Chi2Statistic()
        self.nanChi2.chi2 = np.nan
        self.nanChi2.ndof = 100

        self.maxSteps = 20
        self.name = "testing"
        self.dataName = "fake"
        self.whatToFit = ""  # unneeded, since we're mocking the fitter

        # Mock a Butler so the refObjLoaders have something to call `get()` on.
        self.butler = unittest.mock.Mock(spec=lsst.daf.persistence.Butler)
        self.butler.get.return_value.indexer = DatasetConfig().indexer

        # Mock the association manager and give it access to the ccd list above.
        self.associations = mock.Mock(spec=lsst.jointcal.Associations)
        self.associations.getCcdImageList.return_value = self.ccdImageList

        # a default config to be modified by individual tests
        self.config = lsst.jointcal.jointcal.JointcalConfig()


class TestJointcalIterateFit(JointcalTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        # Mock the fitter and model, so we can force particular
        # return values/exceptions. Default to "good" return values.
        self.fitter = mock.Mock(spec=lsst.jointcal.PhotometryFit)
        self.fitter.computeChi2.return_value = self.goodChi2
        self.fitter.minimize.return_value = MinimizeResult.Converged
        self.model = mock.Mock(spec=lsst.jointcal.SimpleFluxModel)

        self.jointcal = lsst.jointcal.JointcalTask(config=self.config, butler=self.butler)

    def test_iterateFit_success(self):
        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.goodChi2)
        # Once for the for loop, the second time for the rank update.
        self.assertEqual(self.fitter.minimize.call_count, 2)

    def test_iterateFit_writeChi2Outer(self):
        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          self.maxSteps, self.name, self.whatToFit,
                                          dataName=self.dataName)
        self.assertEqual(chi2, self.goodChi2)
        # Once for the for loop, the second time for the rank update.
        self.assertEqual(self.fitter.minimize.call_count, 2)
        # Default config should not call saveChi2Contributions
        self.fitter.saveChi2Contributions.assert_not_called()

    def test_iterateFit_failed(self):
        self.fitter.minimize.return_value = MinimizeResult.Failed

        with self.assertRaises(RuntimeError):
            self.jointcal._iterate_fit(self.associations, self.fitter,
                                       self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(self.fitter.minimize.call_count, 1)

    def test_iterateFit_badFinalChi2(self):
        log = mock.Mock(spec=lsst.log.Log)
        self.jointcal.log = log
        self.fitter.computeChi2.return_value = self.badChi2

        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.badChi2)
        log.info.assert_called_with("%s %s", "Fit completed", self.badChi2)
        log.error.assert_called_with("Potentially bad fit: High chi-squared/ndof.")

    def test_iterateFit_exceedMaxSteps(self):
        log = mock.Mock(spec=lsst.log.Log)
        self.jointcal.log = log
        self.fitter.minimize.return_value = MinimizeResult.Chi2Increased
        maxSteps = 3

        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.goodChi2)
        self.assertEqual(self.fitter.minimize.call_count, maxSteps)
        log.error.assert_called_with("testing failed to converge after %s steps" % maxSteps)

    def test_invalid_model(self):
        self.model.validate.return_value = False
        with(self.assertRaises(ValueError)):
            self.jointcal._logChi2AndValidate(self.associations, self.fitter, self.model)

    def test_nonfinite_chi2(self):
        self.fitter.computeChi2.return_value = self.nanChi2
        with(self.assertRaises(FloatingPointError)):
            self.jointcal._logChi2AndValidate(self.associations, self.fitter, self.model)

    def test_writeChi2(self):
        filename = "somefile"
        self.jointcal._logChi2AndValidate(self.associations, self.fitter, self.model,
                                          writeChi2Name=filename)
        # logChi2AndValidate prepends `config.debugOutputPath` to the filename
        self.fitter.saveChi2Contributions.assert_called_with("./"+filename+"{type}")


class TestJointcalLoadRefCat(JointcalTestBase, lsst.utils.tests.TestCase):

    def _make_fake_refcat(self):
        """Make a fake reference catalog and the bits necessary to use it."""
        center = lsst.geom.SpherePoint(30, -30, lsst.geom.degrees)
        flux = 10
        radius = 1 * lsst.geom.degrees
        filterName = 'fake'

        fakeRefCat = make_fake_refcat(center, flux, filterName)
        fluxField = getRefFluxField(fakeRefCat.schema, filterName)
        returnStruct = lsst.pipe.base.Struct(refCat=fakeRefCat, fluxField=fluxField)
        refObjLoader = mock.Mock(spec=LoadIndexedReferenceObjectsTask)
        refObjLoader.loadSkyCircle.return_value = returnStruct

        return refObjLoader, center, radius, filterName, fakeRefCat

    def test_load_reference_catalog(self):
        refObjLoader, center, radius, filterName, fakeRefCat = self._make_fake_refcat()

        config = lsst.jointcal.jointcal.JointcalConfig()
        config.astrometryReferenceErr = 0.1  # our test refcats don't have coord errors
        jointcal = lsst.jointcal.JointcalTask(config=config, butler=self.butler)

        refCat, fluxField = jointcal._load_reference_catalog(refObjLoader,
                                                             jointcal.astrometryReferenceSelector,
                                                             center,
                                                             radius,
                                                             filterName)
        # operator== isn't implemented for Catalogs, so we have to check like
        # this, in case the records are copied during load.
        self.assertEqual(len(refCat), len(fakeRefCat))
        for r1, r2 in zip(refCat, fakeRefCat):
            self.assertEqual(r1, r2)

    def test_load_reference_catalog_subselect(self):
        """Test that we can select out the one source in the fake refcat
        with a ridiculous S/N cut.
        """
        refObjLoader, center, radius, filterName, fakeRefCat = self._make_fake_refcat()

        config = lsst.jointcal.jointcal.JointcalConfig()
        config.astrometryReferenceErr = 0.1  # our test refcats don't have coord errors
        config.astrometryReferenceSelector.doSignalToNoise = True
        config.astrometryReferenceSelector.signalToNoise.minimum = 1e10
        config.astrometryReferenceSelector.signalToNoise.fluxField = "fake_flux"
        config.astrometryReferenceSelector.signalToNoise.errField = "fake_fluxErr"
        jointcal = lsst.jointcal.JointcalTask(config=config, butler=self.butler)

        refCat, fluxField = jointcal._load_reference_catalog(refObjLoader,
                                                             jointcal.astrometryReferenceSelector,
                                                             center,
                                                             radius,
                                                             filterName)
        self.assertEqual(len(refCat), 0)


class TestJointcalFitModel(JointcalTestBase, lsst.utils.tests.TestCase):
    def test_fit_photometry_writeChi2(self):
        """Test that we are calling saveChi2 with appropriate file prefixes."""
        self.config.photometryModel = "constrainedFlux"
        self.config.writeChi2FilesOuterLoop = True
        jointcal = lsst.jointcal.JointcalTask(config=self.config, butler=self.butler)
        jointcal.focalPlaneBBox = lsst.geom.Box2D()

        # Mock the fitter, so we can pretend it found a good fit
        with mock.patch("lsst.jointcal.PhotometryFit", autospect=True) as fitPatch:
            fitPatch.return_value.computeChi2.return_value = self.goodChi2
            fitPatch.return_value.minimize.return_value = MinimizeResult.Converged

            # config.debugOutputPath is prepended to the filenames that go into saveChi2Contributions
            expected = ["./photometry_init-ModelVisit_chi2", "./photometry_init-Model_chi2",
                        "./photometry_init-Fluxes_chi2", "./photometry_init-ModelFluxes_chi2"]
            expected = [mock.call(x+"-fake{type}") for x in expected]
            jointcal._fit_photometry(self.associations, dataName=self.dataName)
            fitPatch.return_value.saveChi2Contributions.assert_has_calls(expected)

    def test_fit_astrometry_writeChi2(self):
        """Test that we are calling saveChi2 with appropriate file prefixes."""
        self.config.astrometryModel = "constrained"
        self.config.writeChi2FilesOuterLoop = True
        jointcal = lsst.jointcal.JointcalTask(config=self.config, butler=self.butler)
        jointcal.focalPlaneBBox = lsst.geom.Box2D()

        # Mock the fitter, so we can pretend it found a good fit
        fitPatch = mock.patch("lsst.jointcal.AstrometryFit")
        # Mock the projection handler so we don't segfault due to not-fully initialized ccdImages
        projectorPatch = mock.patch("lsst.jointcal.OneTPPerVisitHandler")
        with fitPatch as fit, projectorPatch as projector:
            fit.return_value.computeChi2.return_value = self.goodChi2
            fit.return_value.minimize.return_value = MinimizeResult.Converged
            # return a real ProjectionHandler to keep ConstrainedAstrometryModel() happy
            projector.return_value = lsst.jointcal.IdentityProjectionHandler()

            # config.debugOutputPath is prepended to the filenames that go into saveChi2Contributions
            expected = ["./astrometry_init-DistortionsVisit_chi2", "./astrometry_init-Distortions_chi2",
                        "./astrometry_init-Positions_chi2", "./astrometry_init-DistortionsPositions_chi2"]
            expected = [mock.call(x+"-fake{type}") for x in expected]
            jointcal._fit_astrometry(self.associations, dataName=self.dataName)
            fit.return_value.saveChi2Contributions.assert_has_calls(expected)


class TestComputeBoundingCircle(lsst.utils.tests.TestCase):
    """Tests of Associations.computeBoundingCircle()"""
    def _checkPointsInCircle(self, points, center, radius):
        """Check that all points are within the (center, radius) circle.

        The test is whether the max(points - center) separation is equal to
        (or slightly less than) radius.
        """
        maxSeparation = 0*lsst.geom.degrees
        for point in points:
            maxSeparation = max(maxSeparation, center.separation(point))
        self.assertAnglesAlmostEqual(maxSeparation, radius, maxDiff=3*lsst.geom.arcseconds)
        self.assertLess(maxSeparation, radius)

    def _testPoints(self, ccdImage1, ccdImage2, skyWcs1, skyWcs2, bbox):
        """Fill an Associations object and test that it computes the correct
        bounding circle for the input data.

        Parameters
        ----------
        ccdImage1, ccdImage2 : `lsst.jointcal.CcdImage`
            The CcdImages to add to the Associations object.
        skyWcs1, skyWcs2 : `lsst.afw.geom.SkyWcs`
            The WCS of each of the above images.
        bbox : `lsst.geom.Box2D`
            The ccd bounding box of both images.
        """
        lsst.log.setLevel('jointcal', lsst.log.DEBUG)
        associations = lsst.jointcal.Associations()
        associations.addCcdImage(ccdImage1)
        associations.addCcdImage(ccdImage2)
        associations.computeCommonTangentPoint()

        circle = associations.computeBoundingCircle()
        center = lsst.geom.SpherePoint(circle.getCenter())
        radius = lsst.geom.Angle(circle.getOpeningAngle().asRadians(), lsst.geom.radians)
        points = [lsst.geom.SpherePoint(skyWcs1.pixelToSky(lsst.geom.Point2D(x)))
                  for x in bbox.getCorners()]
        points.extend([lsst.geom.SpherePoint(skyWcs2.pixelToSky(lsst.geom.Point2D(x)))
                      for x in bbox.getCorners()])
        self._checkPointsInCircle(points, center, radius)

    def testPoints(self):
        """Test for points in an "easy" area, far from RA=0 or the poles."""
        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages()
        self._testPoints(struct.ccdImageList[0], struct.ccdImageList[1],
                         struct.skyWcs[0], struct.skyWcs[1], struct.bbox)

    def testPointsRA0(self):
        """Test for CcdImages crossing RA=0; this demonstrates a fix for
        the bug described in DM-19802.
        """
        # Use the same pixel origins as the cfht_minimal data, but put the sky origin at RA=0
        crpix = lsst.geom.Point2D(931.517869, 2438.572109)
        cd = np.array([[5.19513851e-05, -2.81124812e-07],
                      [-3.25186974e-07, -5.19112119e-05]])
        crval1 = lsst.geom.SpherePoint(0.01, -0.01, lsst.geom.degrees)
        crval2 = lsst.geom.SpherePoint(-0.01, 0.01, lsst.geom.degrees)
        wcs1 = lsst.afw.geom.makeSkyWcs(crpix, crval1, cd)
        wcs2 = lsst.afw.geom.makeSkyWcs(crpix, crval2, cd)

        # Put the visit boresights at the WCS origin, for consistency
        visitInfo1 = lsst.afw.image.VisitInfo(exposureId=30577512,
                                              date=DateTime(65321.1),
                                              boresightRaDec=wcs1.getSkyOrigin())
        visitInfo2 = lsst.afw.image.VisitInfo(exposureId=30621144,
                                              date=DateTime(65322.1),
                                              boresightRaDec=wcs1.getSkyOrigin())

        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(fakeWcses=[wcs1, wcs2],
                                                                fakeVisitInfos=[visitInfo1, visitInfo2])
        self._testPoints(struct.ccdImageList[0], struct.ccdImageList[1],
                         struct.skyWcs[0], struct.skyWcs[1], struct.bbox)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
