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

#include <iostream>
#include <cmath>
#include <list>
#include <memory>
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

#include "lsst/log/Log.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Histo2d.h"
#include "lsst/jointcal/Histo4d.h"
#include "lsst/jointcal/FastFinder.h"
#include "lsst/jointcal/ListMatch.h"

namespace {
LOG_LOGGER _log = LOG_GET("lsst.jointcal.ListMatch");
}

namespace lsst {
namespace jointcal {

// cuts.. limits, etc for combinatorial match

/* a Segment is a pair of stars form the same image. it is used for matching starlists */

struct Segment {
    /* data */
    double r, dx, dy;
    std::shared_ptr<const BaseStar> s1, s2;
    int s1rank;

    /* constructor (could set last argument to identity by default)  */
    Segment(std::shared_ptr<const BaseStar> star1, std::shared_ptr<const BaseStar> star2, const int star1Rank,
            const AstrometryTransform &transform) {
        s1rank = star1Rank;
        s1 = std::move(star1);
        s2 = std::move(star2);
        Point P1 = transform.apply(*star1);
        Point P2 = transform.apply(*star2);
        dx = P2.x - P1.x;
        dy = P2.y - P1.y;
        r = sqrt(dx * dx + dy * dy);
    }

    /* arg(other/(*this)) if considered as complex(dx,dy) */
    double relativeAngle(Segment *other) const {
        return atan2(other->dx * dy - dx * other->dy, dx * other->dx + dy * other->dy);
    }

    friend std::ostream &operator<<(std::ostream &stream, const Segment &segment) {
        stream << " dx " << segment.dx << " dy " << segment.dy << " r " << segment.r << std::endl;
        return stream;
    }
};

class SegmentList : public std::list<Segment> {
public:
    //  SegmentList(const BaseStarList &list, const int nStar);
    SegmentList(const BaseStarList &list, int nStar,
                const AstrometryTransform &transform = AstrometryTransformIdentity());
};

using SegmentIterator = std::list<Segment>::iterator;
using SegmentCIterator = std::list<Segment>::const_iterator;

static bool DecreasingLength(const Segment &first, const Segment &second) { return (first.r > second.r); }

SegmentList::SegmentList(const BaseStarList &list, const int nStars, const AstrometryTransform &transform) {
    BaseStarCIterator siStop;

    /* find the fence */
    siStop = list.begin();
    int limit = std::min(nStars, int(list.size())) - 1;  // -1 because test happens after incrementation
    for (int count = 0; count < limit; count++) ++siStop;

    // iterate on star pairs
    int rank = 0;
    for (auto si1 = list.begin(); si1 != siStop; ++si1, rank++)
        for (auto si2 = siStop; si2 != si1; --si2) {
            push_back(Segment(*si1, *si2, rank, transform));
        }
    sort(DecreasingLength); /* allows a break in loops */
}

//#include <pair>

struct SegmentPair : public std::pair<Segment *, Segment *> {
    SegmentPair(Segment *f, Segment *s) : std::pair<Segment *, Segment *>(f, s){};
};

using SegmentPairList = std::list<SegmentPair>;
using SegmentPairListIterator = SegmentPairList::iterator;
using SegmentPairListCIterator = SegmentPairList::const_iterator;

static std::unique_ptr<StarMatchList> MatchListExtract(const SegmentPairList &pairList, int rank1, int rank2,
                                                       const AstrometryTransform &transform) {
    /* first Select in the segment pairs list the ones which make use of star rank1 in segment1
  and star s2 in segment2 */

    std::unique_ptr<StarMatchList> matchList(new StarMatchList);

    for (SegmentPairListCIterator spi = pairList.begin(); spi != pairList.end(); spi++) {
        const SegmentPair &a_pair = *spi;
        if (a_pair.first->s1rank != rank1 || a_pair.second->s1rank != rank2) continue;
        /* now we store as star matches both ends of segment pairs ,
           but only once the beginning of segments because they all have the same,
           given the selection 3 lines above  */
        if (matchList->size() == 0)
            matchList->push_back(StarMatch(transform.apply(*(a_pair.first->s1)), *(a_pair.second->s1),
                                           a_pair.first->s1, a_pair.second->s1));
        /* always store the match at end */
        matchList->push_back(StarMatch(transform.apply(*(a_pair.first->s2)), *(a_pair.second->s2),
                                       a_pair.first->s2, a_pair.second->s2));
    }
    return matchList;
}

static bool DecreasingQuality(const std::unique_ptr<StarMatchList> &first,
                              const std::unique_ptr<StarMatchList> &second) {
    int idiff = first->size() - second->size();
    if (idiff != 0)
        return (idiff > 0);
    else
        return (first->getDist2() < second->getDist2());
}

/* many matching solutions (StarMatchList) will be compared. Store them in a SolList : */

using SolList = std::list<std::unique_ptr<StarMatchList>>;

/* This one searches a general transformation by histogramming the relative size and orientation
of star pairs ( Segment's) built from the 2 lists */

static std::unique_ptr<StarMatchList> ListMatchupRotShift_Old(BaseStarList &list1, BaseStarList &list2,
                                                              const AstrometryTransform &transform,
                                                              const MatchConditions &conditions) {
    SegmentList sList1(list1, conditions.nStarsList1, transform);
    SegmentList sList2(list2, conditions.nStarsList2, AstrometryTransformIdentity());

    /* choose the binning of the histogram so that
       1: ratio = 1 and rotation angle = n * (pi/2) are bin centers. since
       the angle is computed using atan2, its range is [-pi,pi],
       and the histogram range is [-pi-eps, pi-eps], so
       if (angle>pi- angleOffset) angle -= 2*pi before filling.   */
    int nBinsR = 21;
    int nBinsAngle = 180; /* can be divided by 4 */
    double angleOffset = M_PI / nBinsAngle;
    double minRatio = conditions.minSizeRatio();
    double maxRatio = conditions.maxSizeRatio();
    Histo2d histo(nBinsR, minRatio, maxRatio, nBinsAngle, -M_PI - angleOffset, M_PI - angleOffset);

    SegmentIterator segi1, segi2;
    Segment *seg1, *seg2;
    double ratio, angle;
    for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1) {
        seg1 = &(*segi1);
        if (seg1->r == 0) continue;
        for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2) {
            seg2 = &(*segi2);
            /* if one considers the 2 segments as complex numbers z1 and z2, ratio=mod(z1/z2) and angle =
             * arg(z1/z2) */
            /* I did not put a member function in Segment to compute both because we apply a cut on ratio
               before actually
               computing the angle (which involves a call to atan2 (expensive)) */
            ratio = seg2->r / seg1->r;
            if (ratio > maxRatio) continue;
            if (ratio < minRatio) break; /* use th fact that segment lists are sorted by decresing length */
            angle = seg1->relativeAngle(seg2);
            if (angle > M_PI - angleOffset) angle -= 2. * M_PI;
            histo.fill(ratio, angle);
        }
    }
    double binr, bina;
    histo.binWidth(binr, bina);

    SolList Solutions;
    /* now we want to find in the (r,theta) bins that have the highest counts, the star pair
       (one in l1, one in list2) that contribute to the largest number of segment pairs in this bin :
       so, we histogram a couple of integer that uniquely defines the stars, for the segment pairs
       that contribute to the maximum bin. We choose to histogram the rank of s1 of segment 1
       versus the rank of s1 for segment 2 */

    for (int i = 0; i < conditions.maxTrialCount; ++i) {
        double ratioMax, angleMax;
        double maxContent = histo.maxBin(ratioMax, angleMax);
        histo.fill(ratioMax, angleMax, -maxContent);

        if (conditions.printLevel >= 1)
            LOGLS_DEBUG(_log, " valMax " << maxContent << " ratio " << ratioMax << " angle " << angleMax);

        minRatio = ratioMax - binr / 2;
        maxRatio = ratioMax + binr / 2;
        double minAngle = angleMax - bina / 2;
        double maxAngle = angleMax + bina / 2;
        SegmentPairList pairList;
        Histo2d historank(conditions.nStarsList1, 0., conditions.nStarsList1, conditions.nStarsList2, 0.,
                          conditions.nStarsList2);
        /* reloop on segment pairs to select the ones in this specific bin */

        for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1) {
            seg1 = &(*segi1);
            if (seg1->r == 0) continue;
            for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2) {
                seg2 = &(*segi2);
                ratio = seg2->r / seg1->r;
                if (ratio > maxRatio) continue;
                if (ratio < minRatio)
                    break; /* use the fact that segment lists are sorted by decresing length */
                angle = seg1->relativeAngle(seg2);
                if (angle > M_PI - angleOffset) angle -= 2. * M_PI;
                if (angle < minAngle || angle > maxAngle) continue;
                pairList.push_back(SegmentPair(seg1, seg2)); /* store the match */
                historank.fill(seg1->s1rank + 0.5, seg2->s1rank + 0.5);
            }
        }
        for (int iteration = 0; iteration < conditions.maxTrialCount; iteration++) {
            double dr1, dr2;
            double maxval = historank.maxBin(dr1, dr2);
            /* set this bin to zero so that next iteration will find next maximum */
            historank.fill(dr1, dr2, -maxval);
            auto a_list = MatchListExtract(pairList, int(dr1), int(dr2), AstrometryTransformIdentity());
            a_list->refineTransform(conditions.nSigmas);  // mandatory for the sorting fields to be filled
            Solutions.push_back(std::move(a_list));
        }
    } /* end of loop on (r,theta) bins */
    Solutions.sort(DecreasingQuality);
    std::unique_ptr<StarMatchList> best;
    best.swap(*Solutions.begin());
    /* remove the first one from the list */
    Solutions.pop_front();
    if (conditions.printLevel >= 1) {
        LOGLS_DEBUG(_log, "Best solution " << best->computeResidual() << " npairs " << best->size());
        LOGLS_DEBUG(_log, *(best->getTransform()));
        LOGLS_DEBUG(_log, "Chi2 " << best->getChi2() << ',' << " Number of solutions " << Solutions.size());
    }
    return best;
}

/* this matching routine searches brutally a match between lists in
    the 4 parameter space: size ratio, rotation angle, x and y
    shifts. This is done by histogramming where combinations of four
    objets (2 on each list) fall in this 4 parameter space.

    One trick is that rather than using actual offsets, we histogram
    object indices of the combination:
*/

static std::unique_ptr<StarMatchList> ListMatchupRotShift_New(BaseStarList &list1, BaseStarList &list2,
                                                              const AstrometryTransform &transform,
                                                              const MatchConditions &conditions) {
    if (list1.size() <= 4 || list2.size() <= 4) {
        LOGL_FATAL(_log, "ListMatchupRotShift_New : (at least) one of the lists is too short.");
        return nullptr;
    }

    SegmentList sList1(list1, conditions.nStarsList1, transform);
    SegmentList sList2(list2, conditions.nStarsList2, AstrometryTransformIdentity());

    /* choose the binning of the histogram so that
       1: ratio = 1 and rotation angle = n * (pi/2) are bin centers. since
       the angle is computed using atan2, its range is [-pi,pi],
       and the histogram range is [-pi-eps, pi-eps], so
       if (angle>pi- angleOffset) angle -= 2*pi before filling.   */
    int nBinsR = 21;
    int nBinsAngle = 180; /* can be divided by 4 */
    double angleOffset = M_PI / nBinsAngle;
    double minRatio = conditions.minSizeRatio();
    double maxRatio = conditions.maxSizeRatio();
    SparseHisto4d histo(nBinsR, minRatio, maxRatio, nBinsAngle, -M_PI - angleOffset, M_PI - angleOffset,
                        conditions.nStarsList1, 0., conditions.nStarsList1, conditions.nStarsList2, 0.,
                        conditions.nStarsList2, sList1.size() * sList2.size());

    SegmentIterator segi1, segi2;
    Segment *seg1, *seg2;
    double ratio, angle;

    for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1) {
        seg1 = &(*segi1);
        if (seg1->r == 0) continue;
        for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2) {
            seg2 = &(*segi2);
            /* if one considers the 2 segments as complex numbers z1 and z2, ratio=mod(z1/z2) and angle =
             * arg(z1/z2) */
            /* I did not put a member function in Segment to compute both because we apply a cut on ratio
               before actually
               computing the angle (which involves a call to atan2 (expensive)) */
            ratio = seg2->r / seg1->r;
            if (ratio > maxRatio) continue;
            if (ratio < minRatio) break; /* use th fact that segment lists are sorted by decresing length */
            angle = seg1->relativeAngle(seg2);
            if (angle > M_PI - angleOffset) angle -= 2. * M_PI;
            histo.fill(ratio, angle, seg1->s1rank + 0.5, seg2->s1rank + 0.5);
        }
    }

    SolList Solutions;
    /* now we find the highest bins of the histogram, and recover the original objects.
       This involves actually re-looping on the combinations, but it is much
       faster that the original histogram filling loop, since we only compute
       angle and ratio for Segments that have the right first object
    */

    int oldMaxContent = 0;

    for (int i = 0; i < 4 * conditions.maxTrialCount;
         ++i)  // leave a limit to make avoid (almost)  infinite loops
    {
        double pars[4];
        int maxContent = histo.maxBin(pars);
        if (maxContent == 0) break;
        if (conditions.printLevel >= 1) {
            LOGLS_DEBUG(_log, "ValMax " << maxContent << " ratio " << pars[0] << " angle " << pars[1]);
        }
        histo.zeroBin(pars);
        if (i > 0) { /* the match possibilities come out in a random order when they have the same content.
                        so, we stop investigating guesses when the content goes down AND the requested search
                        depth
                        (maxTrialCount) is reached */
            if (maxContent < oldMaxContent && i >= conditions.maxTrialCount) break;
        }
        oldMaxContent = maxContent;
        /* reloop on segment pairs to select the ones in this specific bin */
        int rank1L1 = int(pars[2]);
        int rank1L2 = int(pars[3]);
        double minAngle, maxAngle;
        histo.binLimits(pars, 0, minRatio, maxRatio);
        histo.binLimits(pars, 1, minAngle, maxAngle);

        std::unique_ptr<StarMatchList> a_list(new StarMatchList);

        for (segi1 = sList1.begin(); segi1 != sList1.end(); ++segi1) {
            seg1 = &(*segi1);
            if (seg1->s1rank != rank1L1) continue;
            if (seg1->r == 0) continue;
            for (segi2 = sList2.begin(); segi2 != sList2.end(); ++segi2) {
                seg2 = &(*segi2);
                if (seg2->s1rank != rank1L2) continue;
                // push in the list the match corresponding to end number 1 of segments
                if (a_list->size() == 0)
                    a_list->push_back(StarMatch(*(seg1->s1), *(seg2->s1), seg1->s1, seg2->s1));
                ratio = seg2->r / seg1->r;
                if (ratio > maxRatio) continue;
                if (ratio < minRatio)
                    break; /* use the fact that segment lists are sorted by decresing length */
                angle = seg1->relativeAngle(seg2);
                if (angle > M_PI - angleOffset) angle -= 2. * M_PI;
                if (angle < minAngle || angle > maxAngle) continue;
                /* here we have 2 segments which have the right
                   - length ratio
                   - relative angle
                   - first objects (objects on the end number 1).
                   The objects on the end number 2 are the actual matches : */
                a_list->push_back(StarMatch(*(seg1->s2), *(seg2->s2), seg1->s2, seg2->s2));
            }
        }

        // a basic check for sanity of the algorithm :

        if (int(a_list->size()) != maxContent + 1) {
            LOGLS_ERROR(_log, "There is an internal inconsistency in ListMatchupRotShift.");
            LOGLS_ERROR(_log, "maxContent  = " << maxContent);
            LOGLS_ERROR(_log, "matches->size() = " << a_list->size());
        }
        a_list->refineTransform(conditions.nSigmas);
        Solutions.push_back(std::move(a_list));
    }

    if (Solutions.size() == 0) {
        LOGLS_ERROR(_log, "Error In ListMatchup : not a single pair match.");
        LOGLS_ERROR(_log, "Probably, the relative scale of lists is not within bounds.");
        LOGLS_ERROR(_log, "min/max ratios: " << minRatio << ' ' << maxRatio);
        return nullptr;
    }

    Solutions.sort(DecreasingQuality);
    std::unique_ptr<StarMatchList> best;
    best.swap(*Solutions.begin());
    /* remove the first one from the list */
    Solutions.pop_front();
    if (conditions.printLevel >= 1) {
        LOGLS_INFO(_log, "Best solution " << best->computeResidual() << " npairs " << best->size());
        LOGLS_INFO(_log, *(best->getTransform()));
        LOGLS_INFO(_log, "Chi2 " << best->getChi2() << ", Number of solutions " << Solutions.size());
    }
    return best;
}

static std::unique_ptr<StarMatchList> ListMatchupRotShift(BaseStarList &list1, BaseStarList &list2,
                                                          const AstrometryTransform &transform,
                                                          const MatchConditions &conditions) {
    if (conditions.algorithm == 1)
        return ListMatchupRotShift_Old(list1, list2, transform, conditions);
    else
        return ListMatchupRotShift_New(list1, list2, transform, conditions);
}

std::unique_ptr<StarMatchList> matchSearchRotShift(BaseStarList &list1, BaseStarList &list2,
                                                   const MatchConditions &conditions) {
    list1.fluxSort();
    list2.fluxSort();

    return ListMatchupRotShift(list1, list2, AstrometryTransformIdentity(), conditions);
}

std::unique_ptr<StarMatchList> matchSearchRotShiftFlip(BaseStarList &list1, BaseStarList &list2,
                                                       const MatchConditions &conditions) {
    list1.fluxSort();
    list2.fluxSort();

    AstrometryTransformLinear flip(0, 0, 1, 0, 0, -1);
    std::unique_ptr<StarMatchList> flipped(ListMatchupRotShift(list1, list2, flip, conditions));
    std::unique_ptr<StarMatchList> unflipped(
            ListMatchupRotShift(list1, list2, AstrometryTransformIdentity(), conditions));
    if (!flipped || !unflipped) return std::unique_ptr<StarMatchList>(nullptr);
    if (conditions.printLevel >= 1) {
        LOGLS_DEBUG(_log,
                    "unflipped Residual " << unflipped->computeResidual() << " nused " << unflipped->size());
        LOGLS_DEBUG(_log, "flipped Residual " << flipped->computeResidual() << " nused " << flipped->size());
    }
    if (DecreasingQuality(flipped, unflipped)) {
        if (conditions.printLevel >= 1) LOGL_DEBUG(_log, "Keeping flipped solution.");
        // One should NOT apply the flip to the result because the matchlist
        // (even the flipped one) contains the actual coordinates of stars.
        // MatchListExtract is always called with AstrometryTransformIdentity() as last parameter
        return flipped;
    } else {
        if (conditions.printLevel >= 1) LOGL_DEBUG(_log, "Keeping unflipped solution.");
        return unflipped;
    }
}

#ifdef STORAGE
// timing : 2.5 s for l1 of 1862 objects  and l2 of 2617 objects
std::unique_ptr<AstrometryTransformLinear> listMatchupShift(const BaseStarList &list1,
                                                            const BaseStarList &list2,
                                                            const AstrometryTransform &transform,
                                                            double maxShift) {
    int ncomb = list1.size() * list2.size();
    if (!ncomb) return nullptr;
    int nx;
    if (ncomb > 10000)
        nx = 100;
    else
        nx = (int)sqrt(ncomb);

    Histo2d histo(nx, -maxShift, maxShift, nx, -maxShift, maxShift);

    BaseStarCIterator s1, s2;
    double x1, y1;
    for (s1 = list1.begin(); s1 != list1.end(); ++s1) {
        transform.apply((*s1)->x, (*s1)->y, x1, y1);
        for (s2 = list2.begin(); s2 != list2.end(); ++s2) {
            histo.fill((*s2)->x - x1, (*s2)->y - y1);
        }
    }
    double dx = 0, dy = 0;
    histo.maxBin(dx, dy);
    return std::unique_ptr<AstrometryTransformLinear>(new AstrometryTransformLinearShift(dx, dy));
}
#endif /*STORAGE*/

// timing : 140 ms for l1 of 1862 objects  and l2 of 2617 objects (450 MHz, "-O4") maxShift = 200.
std::unique_ptr<AstrometryTransformLinear> listMatchupShift(const BaseStarList &list1,
                                                            const BaseStarList &list2,
                                                            const AstrometryTransform &transform,
                                                            double maxShift, double binSize) {
    int nx;
    if (binSize == 0) {
        int ncomb = list1.size() * list2.size();
        if (ncomb > 10000)
            nx = 100;
        else
            nx = (int)sqrt(double(ncomb));
        if (!ncomb) return std::unique_ptr<AstrometryTransformLinear>(nullptr);
    } else
        nx = int(2 * maxShift / binSize + 0.5);

    Histo2d histo(nx, -maxShift, maxShift, nx, -maxShift, maxShift);
    double binSizeNew = 2 * maxShift / nx;

    BaseStarCIterator s1;
    FastFinder finder(list2);
    double x1, y1;
    for (s1 = list1.begin(); s1 != list1.end(); ++s1) {
        transform.apply((*s1)->x, (*s1)->y, x1, y1);
        FastFinder::Iterator it = finder.beginScan(Point(x1, y1), maxShift);
        while (*it) {
            auto s2 = *it;
            histo.fill(s2->x - x1, s2->y - y1);
            ++it;
        }
    }
    SolList Solutions;
    for (int i = 0; i < 4; ++i) {
        double dx = 0, dy = 0;
        double count = histo.maxBin(dx, dy);
        histo.fill(dx, dy, -count);  // zero the maxbin
        AstrometryTransformLinearShift shift(dx, dy);
        auto newGuess = compose(shift, transform);
        auto raw_matches = listMatchCollect(list1, list2, newGuess.get(), binSizeNew);
        std::unique_ptr<StarMatchList> matches(new StarMatchList);
        raw_matches->applyTransform(*matches, &transform);
        matches->setTransformOrder(1);
        matches->refineTransform(3.);
        Solutions.push_back(std::move(matches));
    }
    Solutions.sort(DecreasingQuality);
    std::unique_ptr<AstrometryTransformLinear> best(
            new AstrometryTransformLinear(*std::const_pointer_cast<AstrometryTransformLinear>(
                    std::dynamic_pointer_cast<const AstrometryTransformLinear>(
                            Solutions.front()->getTransform()))));
    return best;
}

#ifdef STORAGE

// this is the old fashioned way...

std::unique_ptr<StarMatchList> listMatchCollect_Slow(const BaseStarList &list1, const BaseStarList &list2,
                                                     const AstrometryTransform *guess, const double maxDist) {
    std::unique_ptr<StarMatchList> matches(new StarMatchList);
    /****** Collect ***********/
    for (BaseStarCIterator si = list1.begin(); si != list1.end(); ++si) {
        const Point *p1 = (*si);
        const Point p2 = guess->apply(*p1);
        const BaseStar *neighbour = list2.findClosest(p2);
        if (!neighbour) continue;
        double distance = p2.Distance(*neighbour);
        if (distance < maxDist) {
            matches->push_back(StarMatch(*p1, *neighbour, *si, neighbour));
            // assign the distance, since we have it in hand:
            matches->back().distance = distance;
        }
    }
    return matches;
}
#endif

// here is the real active routine:

std::unique_ptr<StarMatchList> listMatchCollect(const BaseStarList &list1, const BaseStarList &list2,
                                                const AstrometryTransform *guess, const double maxDist) {
    std::unique_ptr<StarMatchList> matches(new StarMatchList);
    /****** Collect ***********/
    FastFinder finder(list2);
    for (BaseStarCIterator si = list1.begin(); si != list1.end(); ++si) {
        auto p1 = (*si);
        Point p2 = guess->apply(*p1);
        auto neighbour = finder.findClosest(p2, maxDist);
        if (!neighbour) continue;
        double distance = p2.Distance(*neighbour);
        if (distance < maxDist) {
            matches->push_back(StarMatch(*p1, *neighbour, p1, neighbour));
            // assign the distance, since we have it in hand:
            matches->back().distance = distance;
        }
    }
    matches->setTransform(guess);

    return matches;
}

#ifdef STORAGE
// unused
//! iteratively collect and fits, with the same transform kind, until the residual increases
std::unique_ptr<StarMatchList> CollectAndFit(const BaseStarList &list1, const BaseStarList &list2,
                                             const AstrometryTransform *guess, const double maxDist) {
    const AstrometryTransform *bestTransform = guess;
    std::unique_ptr<StarMatchList> prevMatch;
    while (true) {
        auto m = listMatchCollect(list1, list2, bestTransform, maxDist);
        m->setTransform(bestTransform);
        m->refineTransform(3.);
        LOGLS_INFO(_log, "Iterating: resid " << m->computeResidual() << " size " << m->size());
        if (!prevMatch ||
            (prevMatch && m->computeResidual() < prevMatch->computeResidual() * 0.999 && m->Chi2() > 0)) {
            prevMatch.swap(m);
            bestTransform = prevMatch->Transform();
        } else {
            break;
        }
    }
    return prevMatch;
}
#endif

std::unique_ptr<StarMatchList> listMatchCollect(const BaseStarList &list1, const BaseStarList &list2,
                                                const double maxDist) {
    std::unique_ptr<StarMatchList> matches(new StarMatchList);
    FastFinder finder(list2);
    for (BaseStarCIterator si = list1.begin(); si != list1.end(); ++si) {
        auto p1 = (*si);
        auto neighbour = finder.findClosest(*p1, maxDist);
        if (!neighbour) continue;
        double distance = p1->Distance(*neighbour);
        if (distance < maxDist) {
            matches->push_back(StarMatch(*p1, *neighbour, p1, neighbour));
            // assign the distance, since we have it in hand:
            matches->back().distance = distance;
        }
    }

    matches->setTransform(std::make_shared<AstrometryTransformIdentity>());

    return matches;
}

static bool is_transform_ok(const StarMatchList *match, double pixSizeRatio2, const size_t nmin) {
    if ((fabs(fabs(std::dynamic_pointer_cast<const AstrometryTransformLinear>(match->getTransform())
                           ->determinant()) -
              pixSizeRatio2) /
                 pixSizeRatio2 <
         0.2) &&
        (match->size() > nmin))
        return true;
    LOGL_ERROR(_log, "transform is not ok!");
    match->printTransform();
    return false;
}

// utility to check current transform difference
static double transform_diff(const BaseStarList &List, const AstrometryTransform *T1,
                             const AstrometryTransform *T2) {
    double diff2 = 0;
    FatPoint tf1;
    Point tf2;
    int count = 0;
    for (BaseStarCIterator it = List.begin(); it != List.end(); ++it) {
        const BaseStar &s = **it;
        T1->transformPosAndErrors(s, tf1);
        T2->apply(s, tf2);
        double dx = tf1.x - tf2.x;
        double dy = tf1.y - tf2.y;
        diff2 += (tf1.vy * dx * dx + tf1.vx * dy * dy - 2 * tf1.vxy * dx * dy) /
                 (tf1.vx * tf1.vy - tf1.vxy * tf1.vxy);
        count++;
    }
    if (count) return diff2 / double(count);
    return 0;
}

static double median_distance(const StarMatchList *match, const AstrometryTransform *transform) {
    size_t nstars = match->size();
    std::vector<double> resid(nstars);
    std::vector<double>::iterator ir = resid.begin();
    for (auto it = match->begin(); it != match->end(); ++it, ++ir)
        *ir = sqrt(transform->apply(it->point1).computeDist2(it->point2));
    sort(resid.begin(), resid.end());
    return (nstars & 1) ? resid[nstars / 2] : (resid[nstars / 2 - 1] + resid[nstars / 2]) * 0.5;
}

std::unique_ptr<AstrometryTransform> listMatchCombinatorial(const BaseStarList &List1,
                                                            const BaseStarList &List2,
                                                            const MatchConditions &conditions) {
    BaseStarList list1, list2;
    List1.copyTo(list1);
    list1.fluxSort();
    List2.copyTo(list2);
    list2.fluxSort();

    LOGLS_INFO(_log, "listMatchCombinatorial: find match between " << list1.size() << " and " << list2.size()
                                                                   << " stars...");
    auto match = matchSearchRotShiftFlip(list1, list2, conditions);
    double pixSizeRatio2 = std::pow(conditions.sizeRatio, 2);
    size_t nmin =
            std::min(size_t(10), size_t(std::min(List1.size(), List2.size()) * conditions.minMatchRatio));

    std::unique_ptr<AstrometryTransform> transform;
    if (is_transform_ok(match.get(), pixSizeRatio2, nmin))
        transform = match->getTransform()->clone();
    else {
        LOGL_ERROR(_log, "listMatchCombinatorial: direct transform failed, trying reverse");
        match = matchSearchRotShiftFlip(list2, list1, conditions);
        if (is_transform_ok(match.get(), pixSizeRatio2, nmin))
            transform = match->inverseTransform();
        else {
            LOGL_FATAL(_log, "FAILED");
        }
    }

    if (transform) {
        LOGL_INFO(_log, "FOUND");
        if (conditions.printLevel >= 1) {
            LOGL_DEBUG(_log, " listMatchCombinatorial: found the following transform.");
            LOGLS_DEBUG(_log, *transform);
        }
    } else
        LOGL_ERROR(_log, "listMatchCombinatorial: failed to find a transform");
    return transform;
}

std::unique_ptr<AstrometryTransform> listMatchRefine(const BaseStarList &List1, const BaseStarList &List2,
                                                     std::unique_ptr<AstrometryTransform> transform,
                                                     const int maxOrder) {
    if (!transform) {
        return std::unique_ptr<AstrometryTransform>(nullptr);
    }

    // some hard-coded constants that could go in a param file
    const double brightDist = 2.;  // distance in pixels in a match
    const double fullDist = 4.;    // distance in pixels in a match between entire lists
    const double nSigmas = 3.;     // k-sigma clipping on residuals
    const size_t nStars = 500;     // max number of bright stars to fit

    int order = 1;
    size_t nstarmin = 3;

    BaseStarList list1, list2;
    List1.copyTo(list1);
    list1.fluxSort();
    list1.cutTail(nStars);
    List2.copyTo(list2);
    list2.fluxSort();
    list2.cutTail(nStars);

    auto fullMatch = listMatchCollect(List1, List2, transform.get(), fullDist);
    auto brightMatch = listMatchCollect(list1, list2, transform.get(), brightDist);
    double curChi2 = computeChi2(*brightMatch, *transform) / brightMatch->size();

    LOGLS_INFO(_log, "listMatchRefine: start: med.resid " << median_distance(fullMatch.get(), transform.get())
                                                          << " #match " << fullMatch->size());

    do {  // loop on transform order on full list of stars
        auto curTransform = brightMatch->getTransform()->clone();
        unsigned iter = 0;
        double transDiff;
        do {  // loop on transform diff only on bright stars
            brightMatch->setTransformOrder(order);
            brightMatch->refineTransform(nSigmas);
            transDiff = transform_diff(list1, brightMatch->getTransform().get(), curTransform.get());
            curTransform = brightMatch->getTransform()->clone();
            brightMatch = listMatchCollect(list1, list2, curTransform.get(), brightDist);
        } while (brightMatch->size() > nstarmin && transDiff > 0.05 && ++iter < 5);

        double prevChi2 = curChi2;
        curChi2 = computeChi2(*brightMatch, *curTransform) / brightMatch->size();

        fullMatch = listMatchCollect(List1, List2, curTransform.get(), fullDist);
        LOGLS_INFO(_log, "listMatchRefine: order " << order << " med.resid "
                                                   << median_distance(fullMatch.get(), curTransform.get())
                                                   << " #match " << fullMatch->size());
        if (((prevChi2 - curChi2) > 0.01 * curChi2) && curChi2 > 0) {
            LOGLS_INFO(_log, " listMatchRefine: order " << order << " was a better guess.");
            transform = brightMatch->getTransform()->clone();
        }
        nstarmin = brightMatch->getTransform()->getNpar();
    } while (++order <= maxOrder);

    return transform;
}
}  // namespace jointcal
}  // namespace lsst
