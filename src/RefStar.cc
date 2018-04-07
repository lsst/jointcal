// -*- C++ -*-
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

BaseStarList &ref2Base(RefStarList &starList) { return (BaseStarList &)starList; }

BaseStarList *ref2Base(RefStarList *starList) { return (BaseStarList *)starList; }

const BaseStarList &ref2Base(const RefStarList &starList) { return (const BaseStarList &)starList; }

const BaseStarList *ref2Base(const RefStarList *starList) { return (BaseStarList *)starList; }
}  // namespace jointcal
}  // namespace lsst
