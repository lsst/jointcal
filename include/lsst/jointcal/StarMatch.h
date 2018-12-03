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

#ifndef LSST_JOINTCAL_STAR_MATCH_H
#define LSST_JOINTCAL_STAR_MATCH_H

#include <iostream>
#include <iterator>   // for std::ostream_iterator
#include <algorithm>  // for swap
#include <string>
#include <list>

#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/BaseStar.h"             // class definition used in inlined functions
#include "lsst/jointcal/AstrometryTransform.h"  // inlined function calls AstrometryTransform::apply()

namespace lsst {
namespace jointcal {

/// \file
/// \brief pairs of points
///
/// Used to handles matches of star std::lists (image/image) or image/catalog
/// to fit geometrical and photometric transformations.
/// a pair of stars, usually belonging to different images.
/// One would normally assume that
///   they are the same object of the sky. The object contains basically two 2d points (called
///   later point1 and point2), and
///   two pointers (unused by the class and its satellites), that enable in the end to trace back
///   the stars in the caller data structures. */

//! A hanger for star associations
class StarMatch {
    friend class StarMatchList;

public:
    /* if one sets that private, then fitting routines will be a nightmare. we could set all the Transform
     * classes friend... */
    FatPoint point1, point2;  //!< 2 points
    //!< the Star pointers (the pointer is in fact generic, pointed data is never used).
    std::shared_ptr<const BaseStar> s1, s2;
    double distance;
    double chi2;

    //! constructor.
    /*! gives 2 points (that contain the geometry), plus pointers to the Star objects
      (which are there for user convenience). */
    StarMatch(const FatPoint &point1, const FatPoint &point2, std::shared_ptr<const BaseStar> star1,
              std::shared_ptr<const BaseStar> star2)
            : point1(point1), point2(point2), s1(std::move(star1)), s2(std::move(star2)), distance(0.){};

    // the next one would require that StarMatch knows BaseStar which is not mandatory for StarMatch to work
    // StarMatch(BaseStar *star1, BaseStar *star2) : point1(*star1), point2(*star2), s1(star1), s2(star2) {};

    //! returns the distance from transform(point1) to point2.
    double computeDistance(const AstrometryTransform &transform) const {
        return point2.Distance(transform.apply(point1));
    };

    //! returns the chi2 (using errors in the FatPoint's)
    double computeChi2(const AstrometryTransform &transform) const;

    //! to be used before sorting on distances.
    void setDistance(const AstrometryTransform &transform) { distance = computeDistance(transform); };
    //! returns the value computed by the above one.
    double getDistance() const { return distance; }

    void swap() {
        std::swap(point1, point2);
        std::swap(s1, s2);
    }

    /* comparison that ensures that after a sort, duplicates are next one another */
    explicit StarMatch(){};

    friend std::ostream &operator<<(std::ostream &stream, const StarMatch &Match);

    ~StarMatch() {}

private:
    //  bool operator <  (const StarMatch & other) const { return (s1 > other.s1) ? s1 > other.s1 : s2 >
    //  other.s2;}

    /* for unique to remove duplicates */
    bool operator==(const StarMatch &other) const { return (s1 == other.s1 && s2 == other.s2); };
    bool operator!=(const StarMatch &other) const { return (s1 != other.s1 || s2 != other.s2); };

    friend bool compareStar1(const StarMatch &one, const StarMatch &two);
    friend bool sameStar1(const StarMatch &one, const StarMatch &two);
    friend bool compareStar2(const StarMatch &one, const StarMatch &two);
    friend bool sameStar2(const StarMatch &one, const StarMatch &two);
};

inline bool compareStar1(const StarMatch &one, const StarMatch &two) {
    return ((one.s1 == two.s1) ? (one.distance < two.distance) : (&(*one.s1) > &(*two.s1)));
}

inline bool sameStar1(const StarMatch &one, const StarMatch &two) { return (one.s1 == two.s1); }

inline bool compareStar2(const StarMatch &one, const StarMatch &two) {
    return ((one.s2 == two.s2) ? (one.distance < two.distance) : (&(*one.s2) > &(*two.s2)));
}

inline bool sameStar2(const StarMatch &one, const StarMatch &two) { return (one.s2 == two.s2); }

/* =================================== StarMatchList
 * ============================================================ */

std::ostream &operator<<(std::ostream &stream, const StarMatch &match);

//#ifdef TO_BE_FIXED
typedef ::std::list<StarMatch>::iterator StarMatchIterator;
typedef ::std::list<StarMatch>::const_iterator StarMatchCIterator;
//#endif

//! A std::list of star matches,
/*! To be used as the argument to AstrometryTransform::fit routines. There is as
well a StarMatch::fit routine which fits a polynomial by default,
although the transform may be user-provided. The
StarMatchList::refineTransform is a convenient tool to reject
outliers. Given two catalogs, one can assemble a StarMatchList using
utilities such as listMatchCollect. StarMatchList's have write
capabilities.  NStarMatchList is a generalization of this 2-match to n-matches.
*/

std::ostream &operator<<(std::ostream &stream, const StarMatchList &starMatchList);

class StarMatchList : public std::list<StarMatch> {
public:
    void refineTransform(double nSigmas);

    //! enables to get a transformed StarMatchList. Only positions are transformed, not attached stars. const
    //! routine: "this" remains unchanged.
    void applyTransform(StarMatchList &transformed, const AstrometryTransform *priorTransform,
                        const AstrometryTransform *posteriorTransform = nullptr) const;

    /* constructor */
    StarMatchList() : _order(0), _chi2(0){};

    //! carries out a fit with outlier rejection

    //! enables to access the fitted transformation.
    std::shared_ptr<const AstrometryTransform> getTransform() const { return _transform; }

    //! access to the sum of squared residuals of the last call to refineTransform.
    double getDist2() const { return _dist2; }

    //! access to the chi2 of the last call to refineTransform.
    double getChi2() const { return _chi2; }

    //! returns the order of the used transform
    int getTransformOrder() const { return _order; }

    //! swaps elements 1 and 2 of each starmatch in std::list.
    void swap();

    //! returns the average 1d Residual (last call to refineTransform)
    double computeResidual() const;

    /*! cleans up the std::list of pairs for pairs that share one of their stars, keeping the closest one.
       The distance is computed using transform. which = 1 (2) removes ambiguities
       on the first (second) term of the match. which=3 does both.*/
    unsigned removeAmbiguities(const AstrometryTransform &transform, int which = 3);

    //! sets a transform between the 2 std::lists and deletes the previous or default one.  No fit.
    void setTransform(const AstrometryTransform *transform) { _transform = transform->clone(); }
    //!
    void setTransform(const AstrometryTransform &transform) { _transform = transform.clone(); }
    void setTransform(std::shared_ptr<AstrometryTransform> transform) { _transform = std::move(transform); }

    //! set transform according to the given order.
    void setTransformOrder(int order);

    /*! returns the inverse transform (swap, fit(refineTransform) , and swap).
           The caller should delete the returned pointer. */
    std::unique_ptr<AstrometryTransform> inverseTransform();

    //! Sets the distance (residual) field of all std::list elements. Mandatory before sorting on distances
    void setDistance(const AstrometryTransform &transform);

    //! deletes the tail of the match std::list
    void cutTail(int nKeep);

    //! count the number of elements for which distance is < mindist
    int recoveredNumber(double mindist) const;

    //! print the matching transformation quality (transform, chi2, residual)
    void dumpTransform(std::ostream &stream = std::cout) const;

    ~StarMatchList(){/* should delete the transform.... or use counted refs*/};

private:
    StarMatchList(const StarMatchList &);  // copies nor properly handled
    void operator=(const StarMatchList &);

    int _order;
    double _chi2;
    double _dist2;
    std::shared_ptr<AstrometryTransform> _transform;
};

//! sum of distance squared
double computeDist2(const StarMatchList &S, const AstrometryTransform &transform);

//! the actual chi2
double computeChi2(const StarMatchList &L, const AstrometryTransform &transform);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_STAR_MATCH_H
