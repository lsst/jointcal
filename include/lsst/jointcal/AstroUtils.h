// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_ASTRO_UTILS_H
#define LSST_JOINTCAL_ASTRO_UTILS_H

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/Frame.h"
namespace lsst {
namespace jointcal {

class Frame;
class Gtransfo;
//! Transform a Frame through a Transfo.
Frame applyTransfo(const Frame& inputframe, const Gtransfo& gtransfo, const WhichTransformed which);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTRO_UTILS_H
