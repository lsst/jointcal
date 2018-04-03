// -*- C++ -*-
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

BaseStarList &Ref2Base(RefStarList &starList) { return (BaseStarList &)starList; }

BaseStarList *Ref2Base(RefStarList *starList) { return (BaseStarList *)starList; }

const BaseStarList &Ref2Base(const RefStarList &starList) { return (const BaseStarList &)starList; }

const BaseStarList *Ref2Base(const RefStarList *starList) { return (BaseStarList *)starList; }
}  // namespace jointcal
}  // namespace lsst
