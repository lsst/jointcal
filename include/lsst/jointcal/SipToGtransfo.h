// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIP_TO_GTRANSFO_H
#define LSST_JOINTCAL_SIP_TO_GTRANSFO_H

#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/jointcal/Gtransfo.h"

namespace jointcal = lsst::jointcal;

namespace lsst {
namespace jointcal {

class Frame;

//! Transform the other way around
std::shared_ptr<lsst::afw::geom::SkyWcs> gtransfoToTanWcs(const lsst::jointcal::TanSipPix2RaDec& wcsTransfo,
                                                          const lsst::jointcal::Frame& ccdFrame,
                                                          const bool noLowOrderSipTerms = false);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIP_TO_GTRANSFO_H
