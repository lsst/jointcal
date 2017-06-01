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
    double mag;
    double emag;
    double col;
    int gen;
    double wmag;
    unsigned indexInMatrix;
    int measurementCount;
    const RefStar* _refStar;

    double flux2;
    double fluxErr;
    double fluxErr2;

public:
    FittedStar()
            : BaseStar(),
              mag(-1),
              emag(-1),
              col(0.),
              gen(-1),
              wmag(0),
              indexInMatrix(-1),
              measurementCount(0),
              _refStar(nullptr),
              fluxErr(-1),
              fluxErr2(-1) {}

    FittedStar(const BaseStar& B)
            : BaseStar(B),
              mag(-1),
              emag(-1),
              col(0.),
              gen(-1),
              wmag(0),
              indexInMatrix(0),
              measurementCount(0),
              _refStar(nullptr),
              flux2(-1),
              fluxErr(-1),
              fluxErr2(-1) {}

    //!
    FittedStar(const MeasuredStar& M);

    //!
    void clearBeforeAssoc() {
        indexInMatrix = -1;
        measurementCount = 0;
        _refStar = nullptr;
        wmag = 0;
    }

    //!
    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " mcount: " << measurementCount;
    }

    //!
    int MeasurementCount() const { return measurementCount; }

    //!
    int& MeasurementCount() { return measurementCount; }

    //! derived using available zero points in input images. In the absence ofZP, ZP= 0.
    double Mag() const { return mag; }
    double& Mag() { return mag; }
    double EMag() const { return emag; }
    double& EMag() { return emag; }
    double Col() const { return col; }
    double& Col() { return col; }
    int Generation() const { return gen; }
    int& Generation() { return gen; }

    //!
    void SetMag(double Value) { mag = Value; }

    //! this routine will hopefully soon disappear.
    void AddMagMeasurement(double MagValue, double MagWeight);

    //! index is a value that a fit can set and reread....
    void SetIndexInMatrix(const unsigned& Index) { indexInMatrix = Index; };

    //!
    int IndexInMatrix() const { return indexInMatrix; }

    //! Set the astrometric reference star associated with this star.
    void setRefStar(const RefStar* _refStar);

    //! Get the astrometric reference star associated with this star.
    const RefStar* getRefStar() const { return _refStar; };

    //! getters
    double Flux() const { return flux; }
    double& Flux() { return flux; }
    double FluxErr() const { return fluxErr; }
    double& FluxErr() { return fluxErr; }

    double Flux2() const { return flux2; }
    double& Flux2() { return flux2; }
    double FluxErr2() const { return fluxErr2; }
    double& FluxErr2() { return fluxErr2; }
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
