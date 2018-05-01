# See COPYRIGHT file at the top of the source tree.
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.base import Struct

__all__ = ["JointcalCoaddTaskConfig", "JointcalCoaddTask"]


class JointcalCoaddTaskConfig(MakeCoaddTempExpTask.ConfigClass):
    """Config for JointcalCoaddTask
    """
    def setDefaults(self):
        MakeCoaddTempExpTask.ConfigClass.setDefaults(self)
        self.doApplyUberCal = True


class JointcalCoaddTask(MakeCoaddTempExpTask):

    ConfigClass = JointcalCoaddTaskConfig

    def getCalExp(self, dataRef, bgSubtracted):
        """Return one "calexp" calibrated exposure

        Parameters
        ----------
        dataRef
            a sensor-level data reference
        bgSubtracted
            return calexp with background subtracted? If False get the
            calexp's background background model and add it to the calexp.

        Returns
        -------
        afw.image.ExposureF
            The calibrated exposure. If config.doApplyUberCal, jointcal
            calibrations will be applied to the returned exposure using
            applyMosaicResults.
        """
        exposure = dataRef.get("calexp")

        if not bgSubtracted:
            background = dataRef.get("calexpBackground")
            mi = exposure.getMaskedImage()
            mi += background.getImage()
            del mi
        if not self.config.doApplyUberCal:
            return exposure
        # if we are here, it means that we have to apply the improved calibrations coming from jointcal
        self.log.info("doApplyUberCal is set - Using jointcal updated calibrations")
        self.applyJointcalResultsExposure(dataRef, calexp=exposure)
        return exposure

    def applyJointcalResultsExposure(self, dataRef, calexp=None):
        """Update an Exposure with the Wcs, from meas_jointcal
        (Calib and flux sacling will be also used later).
        If None, the calexp will be loaded from the dataRef.  Otherwise it is
        updated in-place.
        """
        if calexp is None:
            calexp = dataRef.get("calexp")

        wcsCont = dataRef.get("jointcal_wcs")
        calexp.setWcs(wcsCont)

        return Struct(exposure=calexp)
