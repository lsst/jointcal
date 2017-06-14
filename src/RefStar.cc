// -*- C++ -*-
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

// void RefStar::assignRefFluxes(std::vector<double> const &refFlux) {
//     _refFlux.clear();
//     copy(refFlux.begin(), refFlux.end(), back_inserter(_refFlux));
// }

BaseStarList &Ref2Base(RefStarList &This) { return (BaseStarList &)This; }

BaseStarList *Ref2Base(RefStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Ref2Base(const RefStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Ref2Base(const RefStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
