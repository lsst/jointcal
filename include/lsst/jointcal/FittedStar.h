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

#ifndef LSST_JOINTCAL_FITTED_STAR_H
#define LSST_JOINTCAL_FITTED_STAR_H

#include <iostream>
#include <fstream>

#include "Eigen/Core"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

class MeasuredStar;
class RefStar;
class AstrometryTransform;

/*! \file */

//! objects whose position is going to be fitted. Coordinates in Common Tangent Plane.

struct PmBlock {
    // proper motion in x and y. Units depend on how you fit them
    double pmx, pmy;
    double epmx, epmy, epmxy;
    double color;  // OK it is unrelated, but associated in practice
    bool mightMove;

    PmBlock() : pmx(0), pmy(0), epmx(0), epmy(0), epmxy(0), color(0), mightMove(false){};
};

/**
 * The objects which have been measured several times.
 *
 * MeasuredStars from different CcdImages that represent the same on-sky object all point to one FittedStar.
 */
class FittedStar : public BaseStar, public PmBlock {
public:
    FittedStar() : BaseStar(), _indexInMatrix(-1), _measurementCount(0), _refStar(nullptr) {}

    FittedStar(const BaseStar& baseStar)
            : BaseStar(baseStar), _indexInMatrix(0), _measurementCount(0), _refStar(nullptr) {}

    //!
    FittedStar(const MeasuredStar& measuredStar);

    /// No move, allow copy constructor: we may copy the fitted StarLists when associating and matching
    /// catalogs, otherwise Stars should be managed by shared_ptr only.
    FittedStar(FittedStar const&) = default;
    FittedStar(FittedStar&&) = delete;
    FittedStar& operator=(FittedStar const&) = delete;
    FittedStar& operator=(FittedStar&&) = delete;

    //!
    void clearBeforeAssoc() {
        _indexInMatrix = -1;
        _measurementCount = 0;
        _refStar = nullptr;
        _flux = 0;
        _fluxErr = 0;
        _mag = 0;
        _magErr = 0;
    }

    //!
    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " mcount: " << _measurementCount;
    }

    //!
    int getMeasurementCount() const { return _measurementCount; }
    int& getMeasurementCount() { return _measurementCount; }

    /// Add a measuredStar on-sky magnitude.
    void addMagMeasurement(double magValue, double magWeight);

    //! index is a value that a fit can set and reread....
    void setIndexInMatrix(Eigen::Index const index) { _indexInMatrix = index; };

    //!
    Eigen::Index getIndexInMatrix() const { return _indexInMatrix; }

    //! Set the astrometric reference star associated with this star.
    void setRefStar(const RefStar* _refStar);

    //! Get the astrometric reference star associated with this star.
    const RefStar* getRefStar() const { return _refStar; };

private:
    Eigen::Index _indexInMatrix;
    int _measurementCount;
    const RefStar* _refStar;
};

/****** FittedStarList */

//! A list of FittedStar s. Such a list is typically constructed by Associations
class FittedStarList : public StarList<FittedStar> {
public:
    bool inTangentPlaneCoordinates;

    //!
    FittedStarList() { inTangentPlaneCoordinates = true; }
};

typedef FittedStarList::const_iterator FittedStarCIterator;
typedef FittedStarList::iterator FittedStarIterator;

BaseStarList& Fitted2Base(FittedStarList& This);
BaseStarList* Fitted2Base(FittedStarList* This);
const BaseStarList& Fitted2Base(const FittedStarList& This);
const BaseStarList* Fitted2Base(const FittedStarList* This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FITTED_STAR_H
