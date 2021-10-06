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

"""Code to convert jointcal's output WCS models to distortion maps that can be
used by afw CameraGeom.
"""
import numpy as np

from lsst.afw import cameraGeom
import lsst.afw.geom
import astshim as ast
import lsst.log
from lsst.geom import SpherePoint, Point2D, radians


class CameraModel:
    """Convert a jointcal `~lsst.afw.geom.SkyWcs` into a distortion model and
    detector positions (TODO) that can be used by `~lsst.afw.cameraGeom`.

    Because this code only operates on the WCS, it is independent of the
    format of the persisted output (e.g. gen2 separate files vs. gen3 bundled
    visits).

    Parameters
    ----------
    wcsList : `list` [`lsst.afw.geom.SkyWcs`]
        The WCS to use to compute the distortion model from, preferably from
        multiple visits on the same tract.
    detectors : `list` [`int`]
        Detector ids that correspond one-to-one with ``wcsList``.
    camera : `lsst.afw.cameraGeom.Camera`
        The camera these WCS were fit for.
    n : `int`
        Number of points to compute the pixel scale at, along the +y axis.
    """
    def __init__(self, wcsList, detectors, camera, n=100):
        self.wcsList = wcsList
        self.camera = camera
        self.detectors = detectors
        self.maxFocalRadius = self.camera.computeMaxFocalPlaneRadius()
        self.n = n
        # the computed radius and pixel scales
        self.fieldAngle = None  # degrees
        self.radialScale = None  # arcsec
        self.tangentialScale = None  # arcsec
        # the computed values for every input wcs
        self.fieldAngles = None
        self.radialScales = None
        self.tangentialScales = None
        self.fieldAngleStd = None
        self.radialScaleStd = None
        self.tangentialScaleStd = None

        self.log = lsst.log.Log.getLogger("jointcal.cameraGeom.CameraModel")

    def computeDistortionModel(self):
        """Calculate the afw cameraGeom distortion model to be included in an
        on-disk camera model.


        PLACEHOLDER: This may be as simple as running `computePixelScale` and
        then doing a numpy polynomial fit to it for the cameraGeom input.
        However, we need to check details of how that distortion model is
        stored in a Camera.
            e.g.: np.polyfit(self.fieldAngle, self.radialScale, poly_degree))
        """
        raise NotImplementedError("not yet!")

    def computePixelScale(self):
        """Compute the radial and tangential pixel scale by averaging over
        multiple jointcal WCS models.

        Also computes the standard deviation and logs any WCS that are
        significant outliers.
        The calculations are stored in the ``fieldAngle[s]``,
        ``radialScale[s]``, and ``tangentialScale[s]`` member variables.
        """
        self.fieldAngles = []
        self.radialScales = []
        self.tangentialScales = []
        for id, wcs in zip(self.detectors, self.wcsList):
            fieldAngle, radial, tangential = self._computeDetectorPixelScale(id, wcs)
            self.fieldAngles.append(fieldAngle)
            self.radialScales.append(radial)
            self.tangentialScales.append(tangential)
        # TODO: For now, don't worry about small differences in the computed
        # field angles vs. their respective radial/tangential scales, just
        # assume all fieldAngle positions are "close enough" and warn if not.
        self.fieldAngle = np.mean(self.fieldAngles, axis=0)
        self.fieldAngleStd = np.std(self.fieldAngles, axis=0)
        if self.fieldAngleStd.max() > 1e-4:
            self.log.warning("Large stddev in computed field angles between visits (max: %s degree).",
                             self.fieldAngleStd.max())
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        self.radialScale = np.mean(self.radialScales, axis=0)
        self.radialScaleStd = np.std(self.radialScales, axis=0)
        if self.radialScaleStd.max() > 1e-4:
            self.log.warning("Large stddev in computed radial scales between visits (max: %s arcsec).",
                             self.radialScaleStd.max())
        self.tangentialScale = np.mean(self.tangentialScales, axis=0)
        self.tangentialScaleStd = np.std(self.tangentialScales, axis=0)
        if self.tangentialScaleStd.max() > 1e-4:
            self.log.warning("Large stddev in computed tangential scales between visits (max: %s arcsec).",
                             self.tangentialScaleStd.max())

    def computeCameraPixelScale(self, detector_id=30):
        """Compute the radial and tangential pixel scales using the distortion
        model supplied with the camera.

        This is designed to be directly comparable with the results of
        `~CameraModel.computePixelScale`.

        Parameters
        ----------
        detector_id: `int`
            Detector identifier for the detector_id to use for the calculation.

        Returns
        -------
        fieldAngle : `numpy.ndarray`
            Field angles in degrees.
        radialScale : `numpy.ndarray`
            Radial direction pixel scales in arcseconds/pixel.
        tangentialScale : `numpy.ndarray`
            Tangential direction pixel scales in arcseconds/pixel.
        """
        # Make a trivial SkyWcs to get a field angle->sky transform from.
        iwcToSkyWcs = lsst.afw.geom.makeSkyWcs(Point2D(0, 0), SpherePoint(0, 0, radians),
                                               lsst.afw.geom.makeCdMatrix(1.0 * radians, 0 * radians, True))
        iwcToSkyMap = iwcToSkyWcs.getFrameDict().getMapping("PIXELS", "SKY")
        skyFrame = iwcToSkyWcs.getFrameDict().getFrame("SKY")

        # Extract the transforms that are defined just on the camera.
        pixSys = self.camera[detector_id].makeCameraSys(cameraGeom.PIXELS)
        pixelsToFocal = self.camera.getTransform(pixSys, cameraGeom.FOCAL_PLANE)
        focalToField = self.camera.getTransform(cameraGeom.FOCAL_PLANE, cameraGeom.FIELD_ANGLE)

        # Build a SkyWcs that combines each of the above components.
        pixelFrame = ast.Frame(2, "Domain=PIXELS")
        focalFrame = ast.Frame(2, "Domain=FOCAL")
        iwcFrame = ast.Frame(2, "Domain=IWC")
        frameDict = ast.FrameDict(pixelFrame)
        frameDict.addFrame("PIXELS", pixelsToFocal.getMapping(), focalFrame)
        frameDict.addFrame("FOCAL", focalToField.getMapping(), iwcFrame)
        frameDict.addFrame("IWC", iwcToSkyMap, skyFrame)
        wcs = lsst.afw.geom.SkyWcs(frameDict)

        return self._computeDetectorPixelScale(detector_id, wcs)

    def _computeDetectorPixelScale(self, detector_id, wcs):
        """Compute pixel scale in radial and tangential directions as a
        function of field angle.

        Parameters
        ----------
        detector_id: `int`
            Detector identifier for the detector of this wcs.
        wcs : `lsst.afw.geom.SkyWcs`
            Full focal-plane model to compute pixel scale on.

        Returns
        -------
        fieldAngle : `numpy.ndarray`
            Field angles in degrees.
        radialScale : `numpy.ndarray`
            Radial direction pixel scales in arcseconds/pixel.
        tangentialScale : `numpy.ndarray`
            Tangential direction pixel scales in arcseconds/pixel.

        Notes
        -----
        Pixel scales are calculated from finite differences only along the +y
        focal plane direction.
        """
        focalToSky = wcs.getFrameDict().getMapping('FOCAL', 'SKY')
        mmPerPixel = self.camera[detector_id].getPixelSize()

        focalToPixels = wcs.getFrameDict().getMapping('FOCAL', 'PIXELS')
        trans = wcs.getTransform()  # Pixels to Sky as Point2d -> SpherePoint
        boresight = trans.applyForward(Point2D(focalToPixels.applyForward([0, 0])))

        rs = np.linspace(0, self.maxFocalRadius, self.n)  # focal plane units
        fieldAngle = np.zeros_like(rs)
        radialScale = np.zeros_like(rs)
        tangentialScale = np.zeros_like(rs)
        for i, r in enumerate(rs):
            # point on the sky at position r along the focal plane +y axis
            sp1 = SpherePoint(*focalToSky.applyForward(Point2D([0, r])), radians)
            # point on the sky one pixel further along the focal plane +y axis
            sp2 = SpherePoint(*focalToSky.applyForward(Point2D([0, r + mmPerPixel.getY()])), radians)
            # point on the sky one pixel off of the focal plane +y axis at r
            sp3 = SpherePoint(*focalToSky.applyForward(Point2D([mmPerPixel.getX(), r])), radians)
            fieldAngle[i] = boresight.separation(sp1).asDegrees()
            radialScale[i] = sp1.separation(sp2).asArcseconds()
            tangentialScale[i] = sp1.separation(sp3).asArcseconds()
        return fieldAngle, radialScale, tangentialScale
