// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIP_TO_GTRANSFO_H
#define LSST_JOINTCAL_SIP_TO_GTRANSFO_H

#include "lsst/afw/image/TanWcs.h"
#include "lsst/jointcal/Gtransfo.h"

namespace jointcal = lsst::jointcal;

namespace lsst {
namespace jointcal {

class Frame;

//! Transform an afw TanWcs into a Gtransfo
TanSipPix2RaDec convertTanWcs(const std::shared_ptr<lsst::afw::image::TanWcs> wcs);

//! Transform the other way around
PTR(lsst::afw::image::TanWcs)
gtransfoToTanWcs(const lsst::jointcal::TanSipPix2RaDec wcsTransfo, const lsst::jointcal::Frame &ccdFrame,
                 const bool noLowOrderSipTerms = false);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIP_TO_GTRANSFO_H
