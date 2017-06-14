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

//! The objects which have been measured several times. The MeasuredStar s measuring the same object in
//! differenr CcdImage s point to the same FittedStar.
class FittedStar : public BaseStar, public PmBlock {
    friend class PhotometryFit;
    friend class PhotometryFit2;

private:
    double _mag;
    double _emag;
    int _gen;
    double _wmag;
    unsigned _indexInMatrix;
    int _measurementCount;
    const RefStar* _refStar;

    double _flux2;
    double _fluxErr;
    double _fluxErr2;

public:
    FittedStar()
            : BaseStar(),
              _mag(-1),
              _emag(-1),
              _gen(-1),
              _wmag(0),
              _indexInMatrix(-1),
              _measurementCount(0),
              _refStar(nullptr),
              _fluxErr(-1),
              _fluxErr2(-1) {}

    FittedStar(const BaseStar& baseStar)
            : BaseStar(baseStar),
              _mag(-1),
              _emag(-1),
              _gen(-1),
              _wmag(0),
              _indexInMatrix(0),
              _measurementCount(0),
              _refStar(nullptr),
              _flux2(-1),
              _fluxErr(-1),
              _fluxErr2(-1) {}

    //!
    FittedStar(const MeasuredStar& measuredStar);

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
    void setIndexInMatrix(const unsigned& index) { _indexInMatrix = index; };

    //!
    int getIndexInMatrix() const { return _indexInMatrix; }

    //! Set the astrometric reference star associated with this star.
    void setRefStar(const RefStar* _refStar);

    //! Get the astrometric reference star associated with this star.
    const RefStar* getRefStar() const { return _refStar; };

    //! getters
    double getFluxErr() const { return _fluxErr; }
    double getFlux2() const { return _flux2; }
    double getFluxErr2() const { return _fluxErr2; }
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
