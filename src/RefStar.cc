// -*- C++ -*-
#include <algorithm>
#include <cassert>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

BaseStarList &Ref2Base(RefStarList &This) { return (BaseStarList &)This; }

BaseStarList *Ref2Base(RefStarList *This) { return reinterpret_cast<BaseStarList *>(This); }

BaseStarList const &Ref2Base(RefStarList const &This) { return (BaseStarList const &)This; }

BaseStarList const *Ref2Base(RefStarList const *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
