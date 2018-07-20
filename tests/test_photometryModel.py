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
    def setUp(self):
        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination each time we create a cfht camera below.
        lsst.afw.image.utils.resetFilters()

        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100)
        self.ccdImageList = struct.ccdImageList
        self.camera = struct.camera
        self.catalogs = struct.catalogs
        self.instFluxKeyName = struct.instFluxKeyName

        self.firstIndex = 0  # for assignIndices

        # Set to True in the subclass constructor to do the PhotoCalib calculations in magnitudes.
        self.useMagnitude = False

    def _toPhotoCalib(self, ccdImage, catalog):
        """Test converting this object to a PhotoCalib."""
        photoCalib = self.model.toPhotoCalib(ccdImage)
        if self.useMagnitude:
            result = photoCalib.instFluxToMagnitude(catalog, self.instFluxKeyName)
        else:
            result = photoCalib.instFluxToMaggies(catalog, self.instFluxKeyName)

        pixToFocal = ccdImage.getDetector().getTransform(lsst.afw.cameraGeom.PIXELS,
                                                         lsst.afw.cameraGeom.FOCAL_PLANE)
        stars = lsst.jointcal.testUtils.getMeasuredStarsFromCatalog(catalog, pixToFocal)

        expects = np.empty(len(stars))
        for i, star in enumerate(stars):
            expects[i] = self.model.transform(ccdImage, star)
        self.assertFloatsAlmostEqual(result[:, 0], expects, rtol=1e-13)
        # NOTE: don't compare transformed errors, as they will be different:
        # photoCalib incorporates the model error, while jointcal computes the
        # full covariance matrix, from which the model error should be derived.

    def test_freezeErrorTransform(self):
        """After calling freezeErrorTransform(), the error transform is unchanged
        by offsetParams().
        """
        ccdImage = self.ccdImageList[0]
        catalog = self.catalogs[0]
        pixToFocal = ccdImage.getDetector().getTransform(lsst.afw.cameraGeom.PIXELS,
                                                         lsst.afw.cameraGeom.FOCAL_PLANE)
        stars = lsst.jointcal.testUtils.getMeasuredStarsFromCatalog(catalog, pixToFocal)
        star0 = stars[0]

        self.model.offsetParams(self.delta)
        t1 = self.model.transform(ccdImage, star0)
        t1Err = self.model.transformError(ccdImage, star0)
        self.model.freezeErrorTransform()
        self.model.offsetParams(self.delta)
        t2 = self.model.transform(ccdImage, star0)
        t2Err = self.model.transformError(ccdImage, star0)

        self.assertFloatsNotEqual(t1, t2)
        self.assertFloatsEqual(t1Err, t2Err)


class SimplePhotometryModelTestCase(PhotometryModelTestBase):
    def test_getNpar(self):
        result = self.model.getNpar(self.ccdImageList[0])
        self.assertEqual(result, 1)
        result = self.model.getNpar(self.ccdImageList[1])
        self.assertEqual(result, 1)

    def testGetTotalParameters(self):
        result = self.model.getTotalParameters()
        self.assertEqual(result, 2)

    def test_toPhotoCalib(self):
        self._toPhotoCalib(self.ccdImageList[0], self.catalogs[0])
        self._toPhotoCalib(self.ccdImageList[1], self.catalogs[1])


class SimpleFluxModelTestCase(SimplePhotometryModelTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.photometryModels.SimpleFluxModel(self.ccdImageList)
        self.model.assignIndices("", self.firstIndex)  # have to call this once to let offsetParams work.
        self.delta = np.arange(len(self.ccdImageList), dtype=float)*-0.2 + 1


class SimpleMagnitudeModelTestCase(SimplePhotometryModelTestCase, lsst.utils.tests.TestCase):
    def setUp(self):
        super().setUp()
        self.model = lsst.jointcal.photometryModels.SimpleMagnitudeModel(self.ccdImageList)
        self.model.assignIndices("", self.firstIndex)  # have to call this once to let offsetParams work.
        self.delta = np.arange(len(self.ccdImageList), dtype=float)*-0.2 + 1
        self.useMagnitude = True


class ConstrainedPhotometryModelTestCase(PhotometryModelTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(ConstrainedPhotometryModelTestCase, self).setUp()
        self.visitOrder = 3
        self.focalPlaneBBox = self.camera.getFpBBox()
        self.model = lsst.jointcal.photometryModels.ConstrainedPhotometryModel(self.ccdImageList,
                                                                               self.focalPlaneBBox,
                                                                               self.visitOrder)
        # have to call this once to let offsetParams work.
        self.model.assignIndices("Model", self.firstIndex)
        # tweak to get more than just a constant field for the second ccdImage
        self.delta = np.arange(20, dtype=float)*-0.2 + 1
        # but keep the first ccdImage constant, to help distinguish test failures.
        self.delta[:10] = 0.0
        self.delta[0] = -5.0
        self.model.offsetParams(self.delta)

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

    def test_toPhotoCalib(self):
        self._toPhotoCalib(self.ccdImageList[0], self.catalogs[0])
        self._toPhotoCalib(self.ccdImageList[1], self.catalogs[1])

    def test_assignIndices(self):
        """Test that the correct number of indices were assigned.
        Does not check that the internal mappings are assigned the correct
        indices.
        """
        # need at least two sensors to distinguish "Model" from "ModelVisit"
        # NOTE: createTwoFakeCcdImages() always uses the same two visitIds,
        # so there will be 2 visits total here.
        struct1 = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100, seed=100, fakeCcdId=12)
        ccdImageList = struct1.ccdImageList
        struct2 = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100, seed=101, fakeCcdId=13)
        ccdImageList.extend(struct2.ccdImageList)
        camera = struct1.camera  # the camera is the same in both structs
        visitOrder = 3
        focalPlaneBBox = camera.getFpBBox()
        model = lsst.jointcal.photometryModels.ConstrainedPhotometryModel(ccdImageList,
                                                                          focalPlaneBBox,
                                                                          visitOrder)

        # one polynomial per visit, plus one fitted scale for the second chip.
        expect = 2 * getNParametersPolynomial(self.visitOrder) + 1
        index = model.assignIndices("Model", self.firstIndex)
        self.assertEqual(index, expect)

        # one polynomial per visit
        expect = 2 * getNParametersPolynomial(self.visitOrder)
        index = model.assignIndices("ModelVisit", self.firstIndex)
        self.assertEqual(index, expect)

        # one fitted chip
        expect = 1
        index = model.assignIndices("ModelChip", self.firstIndex)
        self.assertEqual(index, expect)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
