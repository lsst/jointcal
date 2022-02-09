// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cassert>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

Point RefStar::applyProperMotion(Point star, double timeDeltaYears) const {
    if (!_properMotion) {
        return star;
    } else {
        return _properMotion->apply(star, timeDeltaYears);
    }
}

BaseStarList &Ref2Base(RefStarList &This) { return (BaseStarList &)This; }

BaseStarList *Ref2Base(RefStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Ref2Base(const RefStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Ref2Base(const RefStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
