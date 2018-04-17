import numpy as np
import os

import unittest
import lsst.utils.tests

import lsst.afw.cameraGeom
import lsst.afw.geom
import lsst.afw.table
import lsst.afw.image
import lsst.afw.image.utils
import lsst.daf.persistence
import lsst.jointcal.ccdImage
import lsst.jointcal.photometryModels
import lsst.jointcal.star

from lsst.afw.cameraGeom.testUtils import DetectorWrapper


class PhotometryModelTestBase():
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination each time we create a cfht camera below.
        lsst.afw.image.utils.resetFilters()

        # Load or fake the necessary metadata for each CcdImage (using ccd=12)
        input_dir = os.path.join(self.data_dir, 'cfht_minimal')
        visits = [849375, 850587]
        ccdId = 12
        self.butler = lsst.daf.persistence.Butler(input_dir)
        wcs = self.butler.get("calexp_wcs", visit=visits[0], ccd=ccdId)
        visitInfo = self.butler.get("calexp_visitInfo", visit=visits[0], ccd=ccdId)
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0), lsst.afw.geom.Point2I(2000, 2000))
        photoCalib1 = lsst.afw.image.PhotoCalib(100.0, 1.0)
        photoCalib2 = lsst.afw.image.PhotoCalib(120.0, 5.0)
        filt = "some filter"
        dw = DetectorWrapper(id=ccdId, bbox=bbox)
        pixToFocal = dw.detector.getTransform(lsst.afw.cameraGeom.PIXELS, lsst.afw.cameraGeom.FOCAL_PLANE)

        self.instFlux = 5.0
        self.instFluxErr = 2.0
        baseStar0 = lsst.jointcal.star.BaseStar(0, 0, 1, 2)
        self.star0 = lsst.jointcal.star.MeasuredStar(baseStar0)
        self.star0.setInstFlux(self.instFlux)
        point = lsst.afw.geom.Point2D(self.star0.x, self.star0.y)
        pointFocal = pixToFocal.applyForward(point)
        self.star0.setXFocal(pointFocal.getX())
        self.star0.setYFocal(pointFocal.getY())
        baseStar1 = lsst.jointcal.star.BaseStar(100, 200, 3, 4)
        self.star1 = lsst.jointcal.star.MeasuredStar(baseStar1)
        self.star1.setInstFlux(self.instFlux)
        point = lsst.afw.geom.Point2D(self.star1.x, self.star1.y)
        pointFocal = pixToFocal.applyForward(point)
        self.star1.setXFocal(pointFocal.getX())
        self.star1.setYFocal(pointFocal.getY())

        self.firstIndex = 0  # for assignIndices

        # build a fake minimal catalog
        self.instFluxKeyName = "SomeFlux"
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        # centroid
        lsst.afw.table.Point2DKey.addFields(self.schema, "centroid", "centroid", "pixels")
        self.xErrKey = self.schema.addField("centroid_xSigma", type="F")
        self.yErrKey = self.schema.addField("centroid_ySigma", type="F")
        # shape
        self.shapeKey = lsst.afw.table.QuadrupoleKey.addFields(self.schema, "shape", "",
                                                               lsst.afw.table.CoordinateType.PIXEL)
        # Put the fake sources in the minimal catalog.
        self.instFluxKey = self.schema.addField(
            self.instFluxKeyName+"_flux", type="D", doc="post-ISR instFlux")
        self.instFluxErrKey = self.schema.addField(self.instFluxKeyName+"_fluxSigma", type="D",
                                                   doc="post-ISR instFlux stddev")
        self.maggiesKey = self.schema.addField(self.instFluxKeyName+"_calFlux", type="D", doc="maggies")
        self.maggiesErrKey = self.schema.addField(self.instFluxKeyName+"_calFluxErr", type="D",
                                                  doc="maggies stddev")
        self.magnitudeKey = self.schema.addField(self.instFluxKeyName+"_mag", type="D", doc="magnitude")
        self.magnitudeErrKey = self.schema.addField(self.instFluxKeyName+"_magErr", type="D",
                                                    doc="magnitude stddev")
        self.table = lsst.afw.table.SourceTable.make(self.schema)
        self.table.defineCentroid('centroid')
        self.table.defineShape('shape')
        self.catalog = lsst.afw.table.SourceCatalog(self.table)
        record = self.catalog.addNew()
        record.set('id', 1)
        record.set('centroid_x', self.star0.x)
        record.set('centroid_y', self.star0.y)
        record.set(self.instFluxKeyName+'_flux', self.instFlux)
        record.set(self.instFluxKeyName+'_fluxSigma', self.instFluxErr)
        record = self.catalog.addNew()
        record.set('id', 2)
        record.set('centroid_x', self.star1.x)
        record.set('centroid_y', self.star1.y)
        record.set(self.instFluxKeyName+'_flux', self.instFlux*1e-9)
        record.set(self.instFluxKeyName+'_fluxSigma', self.instFluxErr)

        ccdImage1 = lsst.jointcal.ccdImage.CcdImage(self.catalog, wcs, visitInfo, bbox, filt, photoCalib1,
                                                    dw.detector, visits[0], ccdId, self.instFluxKeyName)
        ccdImage2 = lsst.jointcal.ccdImage.CcdImage(self.catalog, wcs, visitInfo, bbox, filt, photoCalib2,
                                                    dw.detector, visits[1], ccdId, self.instFluxKeyName)
        self.ccdImageList = [ccdImage1, ccdImage2]

        # so we can access parts of the camera later (e.g. focal plane)
        self.camera = self.butler.get('camera', visit=visits[0], ccd=ccdId)

    def _toPhotoCalib(self, ccdImage):
        """Test converting this object to a PhotoCalib."""
        photoCalib = self.model.toPhotoCalib(ccdImage)
        expect = self.model.transform(ccdImage, self.star1, self.instFlux)
        point = lsst.afw.geom.Point2D(self.star1.x, self.star1.y)
        self.assertFloatsAlmostEqual(photoCalib.instFluxToMaggies(self.instFlux, point), expect)

    def test_freezeErrorTransform(self):
        """After calling freezeErrorTransform(), the error transform is unchanged
        by offsetParams().
        """
        self.model.offsetParams(self.delta)
        ccdImage = self.ccdImageList[0]
        t1 = self.model.transform(ccdImage, self.star0, self.instFlux)
        t1Err = self.model.transformError(ccdImage, self.star0, self.instFluxErr)
        self.model.freezeErrorTransform()
        self.model.offsetParams(self.delta)
        t2 = self.model.transform(ccdImage, self.star0, self.instFlux)
        t2Err = self.model.transformError(ccdImage, self.star0, self.instFluxErr)

        self.assertFloatsNotEqual(t1, t2)
        self.assertFloatsEqual(t1Err, t2Err)


class SimplePhotometryModelTestCase(PhotometryModelTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(SimplePhotometryModelTestCase, self).setUp()
        self.model = lsst.jointcal.photometryModels.SimplePhotometryModel(self.ccdImageList)
        self.model.assignIndices("", self.firstIndex)  # have to call this once to let offsetParams work.
        self.delta = np.arange(len(self.ccdImageList), dtype=float)*-0.2 + 1

    def test_getNpar(self):
        result = self.model.getNpar(self.ccdImageList[0])
        self.assertEqual(result, 1)
        result = self.model.getNpar(self.ccdImageList[1])
        self.assertEqual(result, 1)

    def test_toPhotoCalib(self):
        self._toPhotoCalib(self.ccdImageList[0])
        self._toPhotoCalib(self.ccdImageList[1])


class ConstrainedPhotometryModelTestCase(PhotometryModelTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(ConstrainedPhotometryModelTestCase, self).setUp()
        self.visitOrder = 3
        self.focalPlaneBBox = self.camera.getFpBBox()
        self.model = lsst.jointcal.photometryModels.ConstrainedPhotometryModel(self.ccdImageList,
                                                                               self.focalPlaneBBox,
                                                                               self.visitOrder)
        # have to call this once to let offsetParams work.
        self.model.assignIndices("", self.firstIndex)
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
        expect = (4*5)/2
        result = self.model.getNpar(self.ccdImageList[0])
        self.assertEqual(result, expect)
        result = self.model.getNpar(self.ccdImageList[1])
        self.assertEqual(result, expect)

    def test_toPhotoCalib(self):
        self._toPhotoCalib(self.ccdImageList[0])
        self._toPhotoCalib(self.ccdImageList[1])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
