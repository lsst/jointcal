#ifndef SIPTOGTRANSFO__H
#define SIPTOGTRANSFO__H

#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/simastrom/Gtransfo.h"

namespace simAstrom = lsst::meas::simastrom; 

namespace lsst {
namespace meas {
namespace simastrom {


TanSipPix2RaDec ConvertTanWcs(const boost::shared_ptr<lsst::afw::image::TanWcs>  wcs);
    
}}} // end of namespaces

#endif /* SIPTOGTRANSFO__H */
