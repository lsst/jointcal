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

#ifndef LSST_JOINTCAL_REF_STAR_H
#define LSST_JOINTCAL_REF_STAR_H

#include <vector>
#include <fstream>

#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

//! Objects used as position anchors, typically USNO stars. Coordinate system defined by user. The Common
//! Tangent Plane seems a good idea.
class RefStar : public BaseStar {
public:
    RefStar(double xx, double yy, double flux, double fluxErr) : BaseStar(xx, yy, flux, fluxErr) {}

    /// No move or copy: each RefStar is unique, and should be accessed/managed via shared_ptr.
    RefStar(RefStar const&) = delete;
    RefStar(RefStar&&) = delete;
    RefStar& operator=(RefStar const&) = default;
    RefStar& operator=(RefStar&&) = delete;
};

/****** RefStarList ***********/

class Frame;

// typedef StarList<RefStar> RefStarList;
class RefStarList : public StarList<RefStar> {};

typedef RefStarList::const_iterator RefStarCIterator;
typedef RefStarList::iterator RefStarIterator;

BaseStarList& Ref2Base(RefStarList& This);
BaseStarList* Ref2Base(RefStarList* This);
const BaseStarList& Ref2Base(const RefStarList& This);
const BaseStarList* Ref2Base(const RefStarList* This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_REF_STAR_H
