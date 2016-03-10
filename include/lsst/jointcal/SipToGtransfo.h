#ifndef SIPTOGTRANSFO__H
#define SIPTOGTRANSFO__H

#include "lsst/afw/image/TanWcs.h"
#include "lsst/jointcal/Gtransfo.h"

namespace jointcal = lsst::jointcal;

namespace lsst {
namespace jointcal {

  class Frame;

  //! Transform an afw TanWcs into a Gtransfo
TanSipPix2RaDec ConvertTanWcs(const boost::shared_ptr<lsst::afw::image::TanWcs>  wcs);

//! Transform the other way around
 PTR(lsst::afw::image::TanWcs)
   GtransfoToTanWcs(const lsst::jointcal::TanSipPix2RaDec WcsTransfo,
		    const lsst::jointcal::Frame &CcdFrame,
		    const bool NoLowOrderSipTerms=false);
    
}} // end of namespaces

#endif /* SIPTOGTRANSFO__H */
