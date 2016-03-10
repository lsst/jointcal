# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.base import Struct
import lsst.afw.image as afwImage

class JointcalCoaddTask(MakeCoaddTempExpTask):
    def getCalExp(self, dataRef, bgSubtracted):
        """!Return one "calexp" calibrated exposure

        @param[in] dataRef        a sensor-level data reference
        @param[in] bgSubtracted   return calexp with background subtracted? If False get the
                                  calexp's background background model and add it to the calexp.
        @return calibrated exposure

        If config.doApplyUberCal, meas_mosaic calibrations will be applied to
        the returned exposure using applyMosaicResults.
        """
        exposure = dataRef.get("calexp", immediate=True)

        if not bgSubtracted:
            background = dataRef.get("calexpBackground", immediate=True)
            mi = exposure.getMaskedImage()
            mi += background.getImage()
            del mi
        if not self.config.doApplyUberCal:
            return exposure
        # if we are here, it means that we have to apply the improved calibrations coming from meas_jointcal
        self.log.info("doApplyUberCal is set - Using jointcal updated calibrations")
        self.applyjointcalResults(dataRef, calexp=exposure)
        return exposure
    
    def applyjointcalResults(self, dataRef, calexp=None):
        """Deprecated function to apply the results to an exposure
        Deprecated, because the mosaic results can be applied to more than
        one kind of target, so it's worth changing the name to be specific.
        """
        return self.applyjointcalResultsExposure(dataRef, calexp).exposure
        
    def applyjointcalResultsExposure(self, dataRef, calexp=None):
        """Update an Exposure with the Wcs, from meas_jointcal
        (Calib and flux sacling will be also used later).
        If None, the calexp will be loaded from the dataRef.  Otherwise it is
        updated in-place.
        """
        if calexp is None:
            calexp = dataRef.get("calexp", immediate=True)

        wcsCont = dataRef.get("wcs",  immediate=True)
        wcs = afwImage.TanWcs.cast(wcsCont.getWcs())
        calexp.setWcs(wcs)
        
        return Struct(exposure=calexp)
