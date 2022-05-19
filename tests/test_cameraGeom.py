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

import os.path

import lsst.utils

from lsst.afw.geom import SkyWcs
import lsst.afw.cameraGeom
import lsst.jointcal.cameraGeometry
from lsst.obs.subaru import HyperSuprimeCam


class TestCameraGeom(lsst.utils.tests.TestCase):
    """Test jointcal's cameraGeom code on data from a jointcal HSC run.
    """
    def setUp(self):
        path = os.path.join(os.path.dirname(__file__), "data/output/HSC-R/9615/")
        # Two WCS from the same visit, which should produce near identical results.
        # 49 is near the center of the focal plane
        self.detector1 = 49
        # 37 is on the west side of the focal plane
        self.detector2 = 37
        # 103 is on the NW edge of the focal plane and rotated 90 degrees
        self.detector3 = 103
        # 90 is on the NE corner of the focal plane
        self.detector4 = 90
        # a selection of WCSs to compare the calculations on
        self.wcs1 = SkyWcs.readFits(os.path.join(path, f"jointcal_wcs-0028976-{self.detector1:03}.fits"))
        self.wcs2 = SkyWcs.readFits(os.path.join(path, f"jointcal_wcs-0028976-{self.detector2:03}.fits"))
        wcs3 = SkyWcs.readFits(os.path.join(path, f"jointcal_wcs-0023924-{self.detector3:03}.fits"))
        wcs4 = SkyWcs.readFits(os.path.join(path, f"jointcal_wcs-0023900-{self.detector4:03}.fits"))
        self.wcsList = [self.wcs1, self.wcs2, wcs3, wcs4]
        self.detectors = [self.detector1, self.detector2, self.detector3, self.detector4]

        hsc = HyperSuprimeCam()
        self.camera = hsc.getCamera()
        self.maxFocalPlaneRadius = self.camera.computeMaxFocalPlaneRadius()

    def testComputeDetectorPixelScale(self):
        """Test CameraModel._computeDetectorPixelScale by comparing two
        detectors from the same visit: they should have near identical pixel
        scale variation.
        """
        # don't need the wcs or detector lists for this test
        model = lsst.jointcal.cameraGeometry.CameraModel([], [], self.camera)
        angle1, radial1, tangential1 = model._computeDetectorPixelScale(self.detector1,
                                                                        self.wcs1)
        angle2, radial2, tangential2 = model._computeDetectorPixelScale(self.detector2,
                                                                        self.wcs2)
        # The angles really should be identical: this value is independent of distortion.
        self.assertFloatsEqual(angle1, angle2)
        # These should be identical, but we allow for machine-precision errors
        # due to possible complexity in the distortion model.
        self.assertFloatsAlmostEqual(radial1, radial2)
        self.assertFloatsAlmostEqual(tangential1, tangential2)
        # Check that the HSC pixel scale is approximately correct at the boresight.
        self.assertFloatsAlmostEqual(radial1[0], 0.168, rtol=1e-2)
        self.assertFloatsAlmostEqual(tangential1[0], 0.168, rtol=1e-2)

    def testComputeCameraPixelScale(self):
        """Test CameraModel.computeCameraPixelScale by comparing two
        detectors: they should have near identical pixel scale variation.
        """
        # don't need the wcs or detector lists for this test
        model = lsst.jointcal.cameraGeometry.CameraModel([], [], self.camera)
        angle1, radial1, tangential1 = model.computeCameraPixelScale(self.detector1)
        angle2, radial2, tangential2 = model.computeCameraPixelScale(self.detector2)
        # For the camera calculation, there can be machine-precision
        # differences between different detectors, so "almost" here.
        self.assertFloatsAlmostEqual(angle1, angle2)
        self.assertFloatsAlmostEqual(radial1, radial2)
        self.assertFloatsAlmostEqual(tangential1, tangential2)
        # Check that the HSC pixel scale is approximately correct at the boresight.
        self.assertFloatsAlmostEqual(radial1[0], 0.168, rtol=1e-3)
        self.assertFloatsAlmostEqual(tangential1[0], 0.168, rtol=1e-3)

    def testComputePixelScale(self):
        """Test that the pixel scale calculation produces something sensible.
        """
        model = lsst.jointcal.cameraGeometry.CameraModel(self.wcsList, self.detectors, self.camera)
        model.computePixelScale()
        angle2, radial2, tangential2 = model._computeDetectorPixelScale(self.detector2,
                                                                        self.wcs2)
        # TOOD: If we ever go with a more careful approach to the mean
        # fieldAngle question, this may have to change.
        self.assertFloatsAlmostEqual(model.fieldAngle, angle2, rtol=1e-4)

        # spot check the second array elements
        self.assertFloatsEqual(model.fieldAngles[1], angle2)
        self.assertFloatsEqual(model.radialScales[1], radial2)
        self.assertFloatsEqual(model.tangentialScales[1], tangential2)
