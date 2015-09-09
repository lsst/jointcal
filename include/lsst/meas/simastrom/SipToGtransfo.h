#ifndef SIPTOGTRANSFO__H
#define SIPTOGTRANSFO__H

#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/simastrom/Gtransfo.h"

namespace simAstrom = lsst::meas::simastrom; 

namespace lsst {
namespace meas {
namespace simastrom {

  class Frame;

  //! Transform an afw TanWcs into a Gtransfo
TanSipPix2RaDec ConvertTanWcs(const boost::shared_ptr<lsst::afw::image::TanWcs>  wcs);

//! Transform the other way around
 PTR(lsst::afw::image::TanWcs) 
   GtransfoToSip(const lsst::meas::simastrom::TanSipPix2RaDec WcsTransfo, 
		 const lsst::meas::simastrom::Frame &CcdFrame);

    
}}} // end of namespaces

#endif /* SIPTOGTRANSFO__H */
