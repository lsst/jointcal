// -*- C++ -*-
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

RefStar::RefStar(const BaseStar &baseStar) : BaseStar(baseStar), index(0) {}

double RefStar::Flux(int filter) const {
    if (filter < 0 && filter >= 10) return FP_INFINITE;
    return refFlux[filter];
}

void RefStar::AssignRefFluxes(std::vector<double> const &reffluxes) {
    refFlux.clear();
    copy(reffluxes.begin(), reffluxes.end(), back_inserter(refFlux));
}

BaseStarList &Ref2Base(RefStarList &This) { return (BaseStarList &)This; }

BaseStarList *Ref2Base(RefStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Ref2Base(const RefStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Ref2Base(const RefStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
