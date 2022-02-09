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
#include "lsst/jointcal/ProperMotion.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

/**
 * Objects used as position/flux anchors (e.g. Gaia DR2 stars). Coordinate system should match that of the
 * fittedStars these are associated with; typically the common tangent plane.
 *
 * RefStars should have their proper motion and parallax corrections pre-applied, so that they are at
 * the same epoch as is stored in Associations.
 */
class RefStar : public BaseStar {
public:
    RefStar(double xx, double yy, double flux, double fluxErr) : BaseStar(xx, yy, flux, fluxErr) {}

    /// No copy: each RefStar is unique, and should be accessed/managed via shared_ptr.
    RefStar(RefStar const&) = delete;
    RefStar(RefStar&&) = default;
    RefStar& operator=(RefStar const&) = delete;
    RefStar& operator=(RefStar&&) = default;

    // pybind11 cannot handle unique_ptr as arguments, so provide this for python-level testing.
    void setProperMotion(ProperMotion const& properMotion) {
        _properMotion = std::make_unique<ProperMotion>(properMotion);
    }
    void setProperMotion(std::unique_ptr<ProperMotion const>&& properMotion) {
        _properMotion = std::move(properMotion);
    }

    /**
     * Apply proper motion correction to the input star, returning a star with PM-corrected
     * coordinates and coordinate errors.
     *
     * Star is returned unchanged if no proper motion data is available.
     *
     * @param star The star to correct for this proper motion.
     * @param timeDeltaYears The difference in time from the correction epoch to correct for, in
     * years.
     *
     * @return The star with corrected coordinates.
     */
    Point applyProperMotion(Point star, double timeDeltaYears) const;

private:
    // RefStars are already PM corrected to a common epoch: this is to correct the associated FittedStar
    // to each MeasuredStar's epoch. Not all refcats have PM data: this will be nullptr if no PM data is
    // available for any reason.
    std::unique_ptr<ProperMotion const> _properMotion;
};

/****** RefStarList ***********/

// typedef StarList<RefStar> RefStarList;
class RefStarList : public StarList<RefStar> {};

using RefStarCIterator = RefStarList::const_iterator;
using RefStarIterator = RefStarList::iterator;

BaseStarList& Ref2Base(RefStarList& This);
BaseStarList* Ref2Base(RefStarList* This);
const BaseStarList& Ref2Base(const RefStarList& This);
const BaseStarList* Ref2Base(const RefStarList* This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_REF_STAR_H
