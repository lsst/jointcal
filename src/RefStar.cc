// -*- C++ -*-
#include <algorithm>
#include <cassert>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

BaseStarList &Ref2Base(RefStarList &This) { return (BaseStarList &)This; }

BaseStarList *Ref2Base(RefStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Ref2Base(const RefStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Ref2Base(const RefStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
