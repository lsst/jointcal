// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FITTED_STAR_H
#define LSST_JOINTCAL_FITTED_STAR_H

#include <iostream>
#include <fstream>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

class MeasuredStar;
class RefStar;
class Gtransfo;

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
    FittedStar()
            : BaseStar(),
              _mag(-1),
              _gen(-1),
              _wmag(0),
              _indexInMatrix(-1),
              _measurementCount(0),
              _refStar(nullptr) {}

    FittedStar(BaseStar const& baseStar)
            : BaseStar(baseStar),
              _mag(-1),
              _gen(-1),
              _wmag(0),
              _indexInMatrix(0),
              _measurementCount(0),
              _refStar(nullptr) {}

    //!
    FittedStar(MeasuredStar const& measuredStar);

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
        _wmag = 0;
    }

    //!
    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " mcount: " << _measurementCount;
    }

    //!
    int getMeasurementCount() const { return _measurementCount; }
    int& getMeasurementCount() { return _measurementCount; }

    //! derived using available zero points in input images. In the absence ofZP, ZP= 0.
    double getMag() const { return _mag; }
    int getGeneration() const { return _gen; }

    //!
    void setMag(double mag) { _mag = mag; }

    //! this routine will hopefully soon disappear.
    void addMagMeasurement(double magValue, double magWeight);

    //! index is a value that a fit can set and reread....
    void setIndexInMatrix(unsigned const& index) { _indexInMatrix = index; };

    //!
    int getIndexInMatrix() const { return _indexInMatrix; }

    //! Set the astrometric reference star associated with this star.
    void setRefStar(RefStar const* _refStar);

    //! Get the astrometric reference star associated with this star.
    RefStar const* getRefStar() const { return _refStar; };

private:
    double _mag;
    int _gen;
    double _wmag;
    unsigned _indexInMatrix;
    int _measurementCount;
    RefStar const* _refStar;

    double _fluxErr;
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
BaseStarList const& Fitted2Base(FittedStarList const& This);
BaseStarList const* Fitted2Base(FittedStarList const* This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FITTED_STAR_H
