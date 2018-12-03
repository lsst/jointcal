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

import numpy as np

import unittest
import lsst.utils.tests
import lsst.jointcal.testUtils

import lsst.afw.cameraGeom
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.image
import lsst.afw.image.utils
import lsst.daf.persistence
import lsst.jointcal.ccdImage
import lsst.jointcal.photometryModels
import lsst.jointcal.star


def getNParametersPolynomial(order):
    """Number of parameters in a photometry polynomial model is (d+1)(d+2)/2."""
    return (order + 1)*(order + 2)/2


class PhotometryModelTestBase:
    """Have the sublass also derive from ``lsst.utils.tests.TestCase`` to cause
    unittest to use the test_* methods in this class.
    """
    def setUp(self):
        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination each time we create a cfht camera below.
        lsst.afw.image.utils.resetFilters()

        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100)
        self.ccdImageList = struct.ccdImageList
        self.camera = struct.camera
        self.catalogs = struct.catalogs
        self.fluxFieldName = struct.fluxFieldName

        self.stars = []
        for catalog, ccdImage in zip(self.catalogs, self.ccdImageList):
            pixToFocal = ccdImage.getDetector().getTransform(lsst.afw.cameraGeom.PIXELS,
                                                             lsst.afw.cameraGeom.FOCAL_PLANE)
            self.stars.append(lsst.jointcal.testUtils.getMeasuredStarsFromCatalog(catalog, pixToFocal))

        self.fittedStar = lsst.jointcal.star.FittedStar(self.stars[0][0])
        # Make a refStar at this fittedStar position, but with different
        # flux and fluxErr, so that it does interesting things when subtracted.
        self.refStar = lsst.jointcal.star.RefStar(self.fittedStar.x,
                                                  self.fittedStar.y,
                                                  self.fittedStar.flux + 50,
                                                  self.fittedStar.fluxErr * 0.01,
                                                  [], [])

        self.firstIndex = 0  # for assignIndices

        # Set to True in the subclass constructor to do the PhotoCalib calculations in magnitudes.
        self.useMagnitude = False

    def _toPhotoCalib(self, ccdImage, catalog, stars):
        """Test converting this object to a PhotoCalib."""
        photoCalib = self.model.toPhotoCalib(ccdImage)
        if self.useMagnitude:
            result = photoCalib.instFluxToMagnitude(catalog, self.fluxFieldName)
        else:
            result = photoCalib.instFluxToMaggies(catalog, self.fluxFieldName)

        expects = np.empty(len(stars))
        for i, star in enumerate(stars):
            expects[i] = self.model.transform(ccdImage, star)
        self.assertFloatsAlmostEqual(result[:, 0], expects, rtol=2e-13)
        # NOTE: don't compare transformed errors, as they will be different:
        # photoCalib incorporates the model error, while jointcal computes the
        # full covariance matrix, from which the model error should be derived.

    def test_toPhotoCalib(self):
        self._toPhotoCalib(self.ccdImageList[0], self.catalogs[0], self.stars[0])
        self._toPhotoCalib(self.ccdImageList[1], self.catalogs[1], self.stars[1])

    def test_freezeErrorTransform(self):
        """After calling freezeErrorTransform(), the error transform is unchanged
        by offsetParams().
        """
        ccdImage = self.ccdImageList[0]
        star0 = self.stars[0][0]

        self.model.offsetParams(self.delta)
        t1 = self.model.transform(ccdImage, star0)
        t1Err = self.model.transformError(ccdImage, star0)
        self.model.freezeErrorTransform()
        self.model.offsetParams(self.delta)
        t2 = self.model.transform(ccdImage, star0)
        t2Err = self.model.transformError(ccdImage, star0)

        self.assertFloatsNotEqual(t1, t2)
        self.assertFloatsEqual(t1Err, t2Err)


class FluxTestBase:
    """Have the sublass also derive from ``lsst.utils.tests.TestCase`` to cause
    unittest to use the test_* methods in this class.
    """
    def test_offsetFittedStar(self):
        value = self.fittedStar.flux

        self.model.offsetFittedStar(self.fittedStar, 0)
        self.assertEqual(self.fittedStar.flux, value)

        self.model.offsetFittedStar(self.fittedStar, 1)
        self.assertEqual(self.fittedStar.flux, value-1)

    def test_computeRefResidual(self):
        result = self.model.computeRefResidual(self.fittedStar, self.refStar)
        self.assertEqual(result, self.fittedStar.flux - self.refStar.flux)


class MagnitudeTestBase:
    """Have the sublass also derive from ``lsst.utils.tests.TestCase`` to cause
    unittest to use the test_* methods in this class.
    """
    def test_offsetFittedStar(self):
        value = self.fittedStar.mag

        self.model.offsetFittedStar(self.fittedStar, 0)
        self.assertEqual(self.fittedStar.mag, value)

        self.model.offsetFittedStar(self.fittedStar, 1)
        self.assertEqual(self.fittedStar.mag, value-1)

    def test_computeRefResidual(self):
        result = self.model.computeRefResidual(self.fittedStar, self.refStar)
        self.assertEqual(result, self.fittedStar.mag - self.refStar.mag)


class SimplePhotometryModelTestBase(PhotometryModelTestBase):
    """Have the sublass also derive from ``lsst.utils.tests.TestCase`` to cause
    unittest to use the test_* methods in this class.
    """
    def test_getNpar(self):
        result = self.model.getNpar(self.ccdImageList[0])
        self.assertEqual(result, 1)
        result = self.model.getNpar(self.ccdImageList[1])
        self.assertEqual(result, 1)

    def testGetTotalParameters(self):
        result = self.model.getTotalParameters()
        self.assertEqual(result, 2)


class SimpleFluxModelTestCase(SimplePhotometryModelTestBase, FluxTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.photometryModels.SimpleFluxModel(self.ccdImageList)
        self.model.assignIndices("", self.firstIndex)  # have to call this once to let offsetParams work.
        self.delta = np.arange(len(self.ccdImageList), dtype=float)*-0.2 + 1


class SimpleMagnitudeModelTestCase(SimplePhotometryModelTestBase,
                                   MagnitudeTestBase,
                                   lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.photometryModels.SimpleMagnitudeModel(self.ccdImageList)
        self.model.assignIndices("", self.firstIndex)  # have to call this once to let offsetParams work.
        self.delta = np.arange(len(self.ccdImageList), dtype=float)*-0.2 + 1
        self.useMagnitude = True


class ConstrainedPhotometryModelTestCase(PhotometryModelTestBase):
    def setUp(self):
        super().setUp()
        self.visitOrder = 3
        self.focalPlaneBBox = self.camera.getFpBBox()
        # Amount to shift the parameters to get more than just a constant field
        # for the second ccdImage.
        # Reverse the range so that the low order terms are the largest.
        self.delta = (np.arange(20, dtype=float)*-0.2 + 1)[::-1]
        # but keep the first ccdImage constant, to help distinguish test failures.
        self.delta[:10] = 0.0
        self.delta[0] = -5.0

    def _initModel2(self, Model):
        """
        Initialize self.model2 with 2 fake sensor catalogs. Call after setUp().

        Parameters
        ----------
        Model : `PhotometryModel`-type
            The PhotometryModel-derived class to construct.
        """
        # We need at least two sensors to distinguish "Model" from "ModelVisit"
        # in `test_assignIndices()`.
        # createTwoFakeCcdImages() always uses the same two visitIds,
        # so there will be 2 visits total here.
        struct1 = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100, seed=100, fakeCcdId=12,
                                                                 photoCalibMean1=1e-2,
                                                                 photoCalibMean2=1.2e-2)
        self.ccdImageList2 = struct1.ccdImageList
        struct2 = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100, seed=101, fakeCcdId=13,
                                                                 photoCalibMean1=2.0e-2,
                                                                 photoCalibMean2=2.2e-2)
        self.ccdImageList2.extend(struct2.ccdImageList)
        camera = struct1.camera  # the camera is the same in both structs
        focalPlaneBBox = camera.getFpBBox()
        self.model2 = Model(self.ccdImageList2, focalPlaneBBox, self.visitOrder)

    def test_getNpar(self):
        """
        Order 3 => (3+1)*(3+2))/2 = 10 parameters,
        and the chip map is fixed (only one ccd), so does not contribute.
        """
        expect = getNParametersPolynomial(self.visitOrder)
        result = self.model.getNpar(self.ccdImageList[0])
        self.assertEqual(result, expect)
        result = self.model.getNpar(self.ccdImageList[1])
        self.assertEqual(result, expect)

    def testGetTotalParameters(self):
        """Two visits, one (fixed) ccd."""
        expect = getNParametersPolynomial(self.visitOrder) * 2
        result = self.model.getTotalParameters()
        self.assertEqual(result, expect)

    def test_assignIndices(self):
        """Test that the correct number of indices were assigned.
        Does not check that the internal mappings are assigned the correct
        indices.
        """
        # one polynomial per visit, plus one fitted scale for the second chip.
        expect = 2 * getNParametersPolynomial(self.visitOrder) + 1
        index = self.model2.assignIndices("Model", self.firstIndex)
        self.assertEqual(index, expect)

        # one polynomial per visit
        expect = 2 * getNParametersPolynomial(self.visitOrder)
        index = self.model2.assignIndices("ModelVisit", self.firstIndex)
        self.assertEqual(index, expect)

        # one fitted chip
        expect = 1
        index = self.model2.assignIndices("ModelChip", self.firstIndex)
        self.assertEqual(index, expect)

    def _testConstructor(self, expectVisit, expectChips):
        """Post-construction, the ChipTransfos should be the PhotoCalib mean of
        the first visit's ccds, and the VisitTransfos should be the identity.
        """
        # Identify to the model that we're fitting both components.
        self.model2.assignIndices("Model", self.firstIndex)

        # check the visitMappings
        for ccdImage in self.ccdImageList2:
            result = self.model2.getMapping(ccdImage).getVisitMapping().getTransform().getParameters()
            self.assertFloatsEqual(result, expectVisit, msg=ccdImage.getName())

        # check the chipMappings
        for ccdImage, expect in zip(self.ccdImageList2, expectChips):
            result = self.model2.getMapping(ccdImage).getChipMapping().getTransform().getParameters()
            # almost equal because log() may have been involved in the math
            self.assertFloatsAlmostEqual(result, expect, msg=ccdImage.getName())


class ConstrainedFluxModelTestCase(ConstrainedPhotometryModelTestCase,
                                   FluxTestBase,
                                   lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.ConstrainedFluxModel(self.ccdImageList,
                                                        self.focalPlaneBBox,
                                                        self.visitOrder)
        # have to call this once to let offsetParams work.
        self.model.assignIndices("Model", self.firstIndex)
        self.model.offsetParams(self.delta)

        self._initModel2(lsst.jointcal.ConstrainedFluxModel)

    def testConstructor(self):
        expectVisit = np.zeros(int(getNParametersPolynomial(self.visitOrder)))
        expectVisit[0] = 1
        # chipMappings are fixed per-chip, and thus are
        # shared between the first pair and second pair of fake ccdImages
        expectChips = [self.ccdImageList2[0].getPhotoCalib().getCalibrationMean(),
                       self.ccdImageList2[0].getPhotoCalib().getCalibrationMean(),
                       self.ccdImageList2[2].getPhotoCalib().getCalibrationMean(),
                       self.ccdImageList2[2].getPhotoCalib().getCalibrationMean()]
        self._testConstructor(expectVisit, expectChips)


class ConstrainedMagnitudeModelTestCase(ConstrainedPhotometryModelTestCase,
                                        MagnitudeTestBase,
                                        lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.ConstrainedMagnitudeModel(self.ccdImageList,
                                                             self.focalPlaneBBox,
                                                             self.visitOrder)
        # have to call this once to let offsetParams work.
        self.model.assignIndices("Model", self.firstIndex)
        self.model.offsetParams(self.delta)

        self._initModel2(lsst.jointcal.ConstrainedMagnitudeModel)

        self.useMagnitude = True

    def testConstructor(self):
        expectVisit = np.zeros(int(getNParametersPolynomial(self.visitOrder)))

        def fluxToMag(flux):
            return -2.5*np.log10(flux)

        # chipMappings are fixed per-chip, and thus are
        # shared between the first pair and second pair of fake ccdImages
        expectChips = [fluxToMag(self.ccdImageList2[0].getPhotoCalib().getCalibrationMean()),
                       fluxToMag(self.ccdImageList2[0].getPhotoCalib().getCalibrationMean()),
                       fluxToMag(self.ccdImageList2[2].getPhotoCalib().getCalibrationMean()),
                       fluxToMag(self.ccdImageList2[2].getPhotoCalib().getCalibrationMean())]
        self._testConstructor(expectVisit, expectChips)

    def test_checkPositiveOnBBox(self):
        self.assertTrue(self.model.checkPositiveOnBBox(self.ccdImageList[0]))
        self.assertTrue(self.model.checkPositiveOnBBox(self.ccdImageList[1]))

        # make a model that is negative all over
        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100, seed=100, fakeCcdId=12,
                                                                photoCalibMean1=1000,
                                                                photoCalibMean2=1200)
        model = lsst.jointcal.ConstrainedMagnitudeModel(struct.ccdImageList,
                                                        struct.camera.getFpBBox(),
                                                        self.visitOrder)
        self.assertFalse(model.checkPositiveOnBBox(struct.ccdImageList[0]))

    def test_validate(self):
        self.assertTrue(self.model.validate(self.ccdImageList))
        # Make the model go negative
        self.model.offsetParams(-3*self.delta)
        self.assertFalse(self.model.validate(self.ccdImageList))


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
