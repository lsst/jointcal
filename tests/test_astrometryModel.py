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

"""Tests of astrometryModels (simple, constrained).

Includes tests of producing a Wcs from a model.
"""
import itertools
import os
import numpy as np
import shutil

import unittest
import lsst.utils.tests

import lsst.afw.cameraGeom
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.image
import lsst.geom
import lsst.log

import lsst.jointcal
from lsst.jointcal.jointcal import (make_schema_table, get_sourceTable_visit_columns,
                                    extract_detector_catalog_from_visit_catalog)
from jointcalTestBase import importRepository


def getNParametersPolynomial(order):
    """Number of parameters in an astrometry polynomial model is 2 * (d+1)(d+2)/2."""
    return (order + 1)*(order + 2)


class AstrometryModelTestBase:
    """Test the jointcal AstrometryModel concrete classes, using CFHT data.
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.dataDir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")

        refcatPath = os.path.join(cls.dataDir, "cfht")
        refcats = {"gaia_dr2_20200414": os.path.join(refcatPath, "gaia_dr2_20200414.ecsv"),
                   "ps1_pv3_3pi_20170110": os.path.join(refcatPath, "ps1_pv3_3pi_20170110.ecsv"),
                   "sdss_dr9_fink_v5b": os.path.join(refcatPath, "sdss-dr9-fink-v5b.ecsv")}
        # Share one repo, since none of these tests write anything.
        cls.repopath = importRepository("lsst.obs.cfht.MegaPrime",
                                        os.path.join(cls.dataDir, 'cfht/repo/'),
                                        os.path.join(cls.dataDir, 'cfht/exports.yaml'),
                                        refcats=refcats,
                                        refcatPath=refcatPath)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.repopath, ignore_errors=True)

    def setUp(self):
        np.random.seed(200)

        # DEBUG messages can help track down failures.
        logger = lsst.log.Log.getLogger('lsst.jointcal')
        logger.setLevel(lsst.log.DEBUG)

        # Append `msg` arguments to assert failures.
        self.longMessage = True
        # absolute tolerance on positional errors of 10 micro-arcsecond
        self.atol = 10.0 / (60 * 60 * 1e6)

        # Maximum difference (see assertPairsAlmostEqual) for round-trip
        # testing of the inverse for models 1 (simpler) and 2 (more.
        # Replace either one for models that don't have as accurate an inverse.
        self.inverseMaxDiff1 = 1e-5
        self.inverseMaxDiff2 = 1e-5

        self.firstIndex = 0  # for assignIndices
        matchCut = 2.0  # arcseconds
        minMeasurements = 2  # accept all star pairs.

        jointcalControl = lsst.jointcal.JointcalControl("flux")
        self.associations = lsst.jointcal.Associations()
        config = lsst.jointcal.JointcalConfig()
        config.load(os.path.join(os.path.dirname(__file__), "config/config.py"))
        sourceSelector = config.sourceSelector.target(config=config.sourceSelector['science'])

        # jointcal's cfht test data has 6 ccds and 2 visits.
        self.visits = [849375, 850587]
        self.detectors = [12, 13, 14, 21, 22, 23]
        self.badVisit = -12345
        self.badCcd = 888

        butler = lsst.daf.butler.Butler(self.repopath, collections='singleFrame', instrument="MegaPrime")

        self.catalogs = []
        self.ccdImageList = []
        table = make_schema_table()
        inColumns = butler.get("sourceTable_visit", visit=self.visits[0])
        columns, ixxColumns = get_sourceTable_visit_columns(inColumns, config, sourceSelector)
        catalogs = {v: sourceSelector.run(butler.get('sourceTable_visit',
                                                     visit=v,
                                                     parameters={'columns': columns})) for v in self.visits}
        for (visit, detector) in itertools.product(self.visits, self.detectors):
            goodSrc = extract_detector_catalog_from_visit_catalog(table,
                                                                  catalogs[visit].sourceCat,
                                                                  detector,
                                                                  ixxColumns,
                                                                  config.sourceFluxType,
                                                                  logger)
            dataId = {"detector": detector, "visit": visit}
            visitInfo = butler.get('calexp.visitInfo', dataId=dataId)
            detector = butler.get('calexp.detector', dataId=dataId)
            ccdId = detector.getId()
            wcs = butler.get('calexp.wcs', dataId=dataId)
            bbox = butler.get('calexp.bbox', dataId=dataId)
            filt = butler.get('calexp.filter', dataId=dataId)
            filterName = filt.physicalLabel
            photoCalib = lsst.afw.image.PhotoCalib(100.0, 1.0)

            self.catalogs.append(goodSrc)
            self.associations.createCcdImage(goodSrc,
                                             wcs,
                                             visitInfo,
                                             bbox,
                                             filterName,
                                             photoCalib,
                                             detector,
                                             visit,
                                             ccdId,
                                             jointcalControl)

        # Have to set the common tangent point so projectionHandler can use skyToCTP.
        self.associations.computeCommonTangentPoint()

        self.projectionHandler = lsst.jointcal.OneTPPerVisitHandler(self.associations.getCcdImageList())

        self.associations.associateCatalogs(matchCut)
        self.associations.prepareFittedStars(minMeasurements)
        self.associations.deprojectFittedStars()

    def _prepModels(self):
        """Call this after model1 and model2 are created, to call assignIndices,
        and instantiate the fitters.
        """
        posError = 0.02  # in pixels
        # have to call this once or offsetParams will fail because the transform indices aren't defined
        self.model1.assignIndices("Distortions", self.firstIndex)
        self.fitter1 = lsst.jointcal.AstrometryFit(self.associations, self.model1, posError)

        # have to call this once or offsetParams will fail because the transform indices aren't defined
        self.model2.assignIndices("Distortions", self.firstIndex)
        self.fitter2 = lsst.jointcal.AstrometryFit(self.associations, self.model2, posError)

    def testMakeSkyWcsModel1(self):
        self.checkMakeSkyWcsModel(self.model1, self.fitter1, self.inverseMaxDiff1)

    def testMakeSkyWcsModel2(self):
        self.checkMakeSkyWcsModel(self.model2, self.fitter2, self.inverseMaxDiff2)

    def checkMakeSkyWcsModel(self, model, fitter, inverseMaxDiff):
        """Test producing a SkyWcs on a model for every cdImage,
        both post-initialization and after one fitting step.

        Parameters
        ----------
        model : `lsst.jointcal.AstrometryModel`
            The model to test.
        fitter : `lsst.jointcal.FitterBase`
            The fitter to use to step the model to test with new (reasonable) parameters.
        inverseMaxDiff : `float`
            Required accuracy on inverse transform.
            See `lsst.afw.geom.utils.assertPairsAlmostEqual`.

        """
        # first test on as-initialized models
        for ccdImage in self.associations.getCcdImageList():
            self.checkMakeSkyWcsOneCcdImage(model, ccdImage, inverseMaxDiff)

        # now shift the models to non-default, but more reasonable, values by taking one fitting step.
        fitter.minimize("DistortionsVisit")
        fitter.minimize("Distortions")
        for ccdImage in self.associations.getCcdImageList():
            self.checkMakeSkyWcsOneCcdImage(model, ccdImage, inverseMaxDiff)

    def checkMakeSkyWcsOneCcdImage(self, model, ccdImage, inverseMaxDiff):
        """Test converting the model of one ccdImage to a SkyWcs by comparing
        to the original transform at the tangent plane.

        Parameters
        ----------
        model : `lsst.jointcal.AstrometryModel`
            The model to test.
        ccdImage : `lsst.jointcal.CcdImage`
            The ccdImage to extract from the model and test.
        inverseMaxDiff : `float`
            Required accuracy on inverse transform.
            See `lsst.afw.geom.utils.assertPairsAlmostEqual`.
        """
        skyWcs = model.makeSkyWcs(ccdImage)
        skyToTangentPlane = model.getSkyToTangentPlane(ccdImage)
        mapping = model.getMapping(ccdImage)

        bbox = ccdImage.getDetector().getBBox()
        num = 200
        xx = np.linspace(bbox.getMinX(), bbox.getMaxX(), num)
        yy = np.linspace(bbox.getMinY(), bbox.getMaxY(), num)
        points = [lsst.geom.Point2D(*xy) for xy in itertools.product(xx, yy)]

        expects = []
        forwards = []
        inverses = []
        spherePoints = skyWcs.pixelToSky(points)
        inverses = skyWcs.skyToPixel(skyWcs.pixelToSky(points))
        for point, spherePoint in zip(points, spherePoints):
            # TODO: Fix these "Point"s once DM-4044 is done.

            # jointcal's pixel->tangent-plane mapping
            star = lsst.jointcal.BaseStar(point.getX(), point.getY(), 0, 0)
            tpExpect = mapping.transformPosAndErrors(star)
            expects.append(lsst.geom.Point2D(tpExpect.x, tpExpect.y))

            # skywcs takes pixel->sky, and we then have to go sky->tangent-plane
            onSky = lsst.jointcal.BaseStar(spherePoint.getLongitude().asDegrees(),
                                           spherePoint.getLatitude().asDegrees(), 0, 0)
            result = skyToTangentPlane.apply(onSky)
            forwards.append(lsst.geom.Point2D(result.x, result.y))

        self.assertPairListsAlmostEqual(forwards, expects)
        # NOTE: assertPairListsAlmostEqual() compares absolute, not relative,
        # values so the points along the ccd edge may exceed maxDiff while still
        # being "close enough": set `inverseMaxDiff` accordingly.
        self.assertPairListsAlmostEqual(inverses, points, maxDiff=inverseMaxDiff)


class SimpleAstrometryModelTestCase(AstrometryModelTestBase, lsst.utils.tests.TestCase):
    """Test the `SimpleAstrometryModel`, with one mapping per ccd per visit."""
    def setUp(self):
        super().setUp()
        self.order1 = 3
        self.inverseMaxDiff1 = 2e-4
        self.model1 = lsst.jointcal.SimpleAstrometryModel(self.associations.getCcdImageList(),
                                                          self.projectionHandler,
                                                          True,
                                                          order=self.order1)

        self.order2 = 5
        # NOTE: because assertPairListsAlmostEqual tests an absolute
        # difference, we need this to be relatively high to avoid spurious
        # incorrect values.
        # Alternately, further increasing the order of the inverse polynomial
        # in astrometryTransform.toAstMap() can improve the quality of the
        # SkyWcs inverse, but that may not be wise for the more general use
        # case due to the inverse then having too many wiggles.
        self.inverseMaxDiff2 = 2e-2
        self.model2 = lsst.jointcal.SimpleAstrometryModel(self.associations.getCcdImageList(),
                                                          self.projectionHandler,
                                                          False,
                                                          order=self.order2)
        self._prepModels()

    def _testGetNpar(self, model, order):
        for ccdImage in self.associations.getCcdImageList():
            result = model.getNpar(ccdImage)
            self.assertEqual(result, getNParametersPolynomial(order))

    def testGetNpar1(self):
        self._testGetNpar(self.model1, self.order1)

    def testGetNpar2(self):
        self._testGetNpar(self.model2, self.order2)

    def _testGetTotalParameters(self, model, order):
        result = model.getTotalParameters()
        expect = getNParametersPolynomial(order)*len(self.associations.getCcdImageList())
        self.assertEqual(result, expect)

    def testGetTotalParametersModel1(self):
        self._testGetTotalParameters(self.model1, self.order1)

    def testGetTotalParametersModel2(self):
        self._testGetTotalParameters(self.model2, self.order2)


class ConstrainedAstrometryModelTestCase(AstrometryModelTestBase, lsst.utils.tests.TestCase):
    """Test the `ConstrainedAstrometryModel`, with one mapping per ccd and one
    mapping per visit.
    """
    def setUp(self):
        super().setUp()
        self.visitOrder1 = 3
        self.chipOrder1 = 1
        self.inverseMaxDiff1 = 1e-5
        self.model1 = lsst.jointcal.ConstrainedAstrometryModel(self.associations.getCcdImageList(),
                                                               self.projectionHandler,
                                                               chipOrder=self.chipOrder1,
                                                               visitOrder=self.visitOrder1)

        self.visitOrder2 = 5
        self.chipOrder2 = 2
        self.inverseMaxDiff2 = 8e-5
        self.model2 = lsst.jointcal.ConstrainedAstrometryModel(self.associations.getCcdImageList(),
                                                               self.projectionHandler,
                                                               chipOrder=self.chipOrder2,
                                                               visitOrder=self.visitOrder2)
        self._prepModels()

        # 22 is closest to the center of the focal plane in this data, so it is not fit.
        self.fixedCcd = 22

    def _polyParams(self, chipOrder, visitOrder):
        """Number of parameters per polynomial is (d+1)(d+2)/2, summed over
        polynomials, times 2 polynomials per dimension.
        The chip transform is fixed for one chip, so only visitOrder matters
        if chipOrder is None.
        """
        params = getNParametersPolynomial(visitOrder)
        if chipOrder is not None:
            params += getNParametersPolynomial(chipOrder)
        return params

    def _testGetNpar(self, model, chipOrder, visitOrder):
        def checkParams(ccdImage, model, chipOrder, visitOrder):
            result = model.getNpar(ccdImage)
            failMsg = "ccdImage: %s, with chipOrder %s and visitOrder %s"%(ccdImage.getName(),
                                                                           chipOrder,
                                                                           visitOrder)
            self.assertEqual(result, self._polyParams(chipOrder, visitOrder), msg=failMsg)

        for ccdImage in self.associations.getCcdImageList():
            realChipOrder = None if ccdImage.getCcdId() == self.fixedCcd else chipOrder
            checkParams(ccdImage, model, realChipOrder, visitOrder)

    def testGetNpar1(self):
        self._testGetNpar(self.model1, self.chipOrder1, self.visitOrder1)

    def testGetNpar2(self):
        self._testGetNpar(self.model2, self.chipOrder2, self.visitOrder2)

    def _testGetTotalParameters(self, model, chipOrder, visitOrder):
        result = model.getTotalParameters()
        # one sensor is held fixed, hence len(ccds)-1
        expect = getNParametersPolynomial(chipOrder)*(len(self.detectors) - 1) + \
            getNParametersPolynomial(visitOrder)*len(self.visits)
        self.assertEqual(result, expect)

    def testGetTotalParametersModel1(self):
        self._testGetTotalParameters(self.model1, self.chipOrder1, self.visitOrder1)

    def testGetTotalParametersModel2(self):
        self._testGetTotalParameters(self.model2, self.chipOrder2, self.visitOrder2)

    def checkGetChipTransform(self, model):
        # Check valid ccds
        for ccd in self.detectors:
            try:
                model.getChipTransform(ccd)
            except lsst.pex.exceptions.wrappers.InvalidParameterError:
                self.fail("model: {} raised on ccd: {}, but should not have.".format(model, ccd))

        # Check an invalid ccd
        with self.assertRaises(lsst.pex.exceptions.wrappers.InvalidParameterError) as cm:
            model.getChipTransform(self.badCcd)
        errMsg = f"No such chipId: {self.badCcd} among [{', '.join(str(d) for d in self.detectors)}]"
        self.assertIn(errMsg, str(cm.exception))

    def testGetChipTransform(self):
        """getChipTransform should get each known transform, and raise with an
        appropriate message otherwise.
        """
        self.checkGetChipTransform(self.model1)
        self.checkGetChipTransform(self.model2)

    def checkGetVisitTransform(self, model):
        # Check valid visits
        for visit in self.visits:
            try:
                model.getVisitTransform(visit)
            except lsst.pex.exceptions.wrappers.InvalidParameterError:
                self.fail("model: {} raised on visit: {}, but should not have.".format(model, visit))

        # Check an invalid visit
        with self.assertRaises(lsst.pex.exceptions.wrappers.InvalidParameterError) as cm:
            model.getVisitTransform(self.badVisit)
        errMsg = "No such visitId: {} among [{}]".format(self.badVisit,
                                                         ", ".join(str(v) for v in self.visits))
        self.assertIn(errMsg, str(cm.exception))

    def testGetVisitTransform(self):
        """getVisitTransform should get each known transform, and raise with an
        appropriate message otherwise.
        """
        self.checkGetVisitTransform(self.model1)
        self.checkGetVisitTransform(self.model2)

    def testValidate(self):
        """Test that invalid models fail validate(), and that valid ones pass.
        """
        # We need at least 0 degrees of freedom (data - parameters) for the model to be valid.
        # Note: model1 has 70 total parameters (2 visits*20 params + (6-1) sensors*5 params)
        self.assertTrue(self.model1.validate(self.ccdImageList, 0))
        self.assertFalse(self.model1.validate(self.ccdImageList, -1))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
