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

#ifndef LSST_JOINTCAL_LIST_MATCH_H
#define LSST_JOINTCAL_LIST_MATCH_H

#include <string>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarMatch.h"

namespace lsst {
namespace jointcal {

class AstrometryTransform;
class AstrometryTransformLinear;

//! Parameters to be provided to combinatorial searches
struct MatchConditions {
    int nStarsList1, nStarsList2;
    int maxTrialCount;
    double nSigmas;
    double maxShiftX, maxShiftY;
    double sizeRatio, deltaSizeRatio, minMatchRatio;
    int printLevel;
    int algorithm;

    MatchConditions()
            : nStarsList1(70),
              nStarsList2(70),
              maxTrialCount(4),
              nSigmas(3.),
              maxShiftX(50),
              maxShiftY(50),
              sizeRatio(1),
              deltaSizeRatio(0.1 * sizeRatio),
              minMatchRatio(1. / 3.),
              printLevel(0),
              algorithm(2) {}

    double minSizeRatio() const { return sizeRatio - deltaSizeRatio; }
    double maxSizeRatio() const { return sizeRatio + deltaSizeRatio; }
};

/*! \file
    \brief Combinatorial searches for linear transformations to go from
           list1 to list2.

    The following routines search a geometrical transformation that make
two lists of stars to match geometrically as well as possible. They are used
either to match two images of the same sky area, or an image with a catalogue.
They assume that fluxes assigned to stars are actual fluxes, i.e. the brighter
the star, the higher the flux. They however only rely on flux ordering,
not values.
 */

//! searches a geometrical transformation that goes from list1 to list2.
/*!  The found transformation is a field of the returned object, as well as the star pairs
(the matches) that were constructed.  (see StarMatchList class definition for more details).
The various cuts are contained in conditions (see listmatch.h) for its contents.
This routine searches a transformation that involves a shift and a rotation. */

std::unique_ptr<StarMatchList> matchSearchRotShift(BaseStarList &list1, BaseStarList &list2,
                                                   const MatchConditions &conditions);

//! same as above but searches also a flipped solution.

std::unique_ptr<StarMatchList> matchSearchRotShiftFlip(BaseStarList &list1, BaseStarList &list2,
                                                       const MatchConditions &conditions);

//! assembles star matches.
/*! It picks stars in list1, transforms them through guess, and collects
closest star in list2, and builds a match if closer than maxDist). */

std::unique_ptr<StarMatchList> listMatchCollect(const BaseStarList &list1, const BaseStarList &list2,
                                                const AstrometryTransform *guess, double maxDist);

//! same as before except that the transform is the identity

std::unique_ptr<StarMatchList> listMatchCollect(const BaseStarList &list1, const BaseStarList &list2,
                                                double maxDist);

//! searches for a 2 dimensional shift using a very crude histogram method.

std::unique_ptr<AstrometryTransformLinear> listMatchupShift(const BaseStarList &list1,
                                                            const BaseStarList &list2,
                                                            const AstrometryTransform &transform,
                                                            double maxShift, double binSize = 0);

std::unique_ptr<AstrometryTransform> listMatchCombinatorial(
        const BaseStarList &list1, const BaseStarList &list2,
        const MatchConditions &conditions = MatchConditions());
std::unique_ptr<AstrometryTransform> listMatchRefine(const BaseStarList &list1, const BaseStarList &list2,
                                                     std::unique_ptr<AstrometryTransform> transform,
                                                     int maxOrder = 3);

#ifdef DO_WE_NEED_THAT
inline AstrometryTransform *ListMatch(const BaseStarList &list1, const BaseStarList &list2,
                                      const int maxOrder = 3) {
    AstrometryTransform *transform = listMatchCombinatorial(list1, list2);
    transform = listMatchRefine(list1, list2, transform, maxOrder);
    return transform;
}
#endif /*  DO_WE_NEED_THAT */
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_LIST_MATCH_H
