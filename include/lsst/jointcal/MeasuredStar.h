// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_MEASURED_STAR_H
#define LSST_JOINTCAL_MEASURED_STAR_H

#include <iostream>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

class CcdImage;

//! objects measured on actual images. Coordinates and uncertainties are expressed in pixel image frame. Flux
//! expressed in ADU/s.
class MeasuredStar : public BaseStar {
public:
    double mag;
    double wmag;
    double chi2;

private:
    // on-chip flux, in ADU
    double _instFlux;
    double _instFluxErr;

    const CcdImage *_ccdImage;
    std::shared_ptr<const FittedStar> _fittedStar;
    bool _valid;

public:
    //!
    MeasuredStar()
            : BaseStar(), mag(0.), wmag(0.), _instFlux(0.), _instFluxErr(0.), _ccdImage(0), _valid(true) {}

    MeasuredStar(const BaseStar &baseStar)
            : BaseStar(baseStar), mag(0.), wmag(0.), _instFluxErr(0.), _ccdImage(0), _valid(true) {}

    void setFittedStar(std::shared_ptr<FittedStar> fittedStar) {
        if (fittedStar) fittedStar->getMeasurementCount()++;
        _fittedStar = std::move(fittedStar);
    }

    void dump(std::ostream &stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " ccdImage: " << _ccdImage << " valid: " << _valid;
    }

    void setInstFlux(double instFlux) { _instFlux = instFlux; }
    void setInstFluxErr(double instFluxErr) { _instFluxErr = instFluxErr; }

    double getInstFlux() const { return _instFlux; }
    double getInstFluxErr() const { return _instFluxErr; }
    double getMag() const { return mag; }

    //! the inverse of the mag variance
    double getMagWeight() const { return (_instFlux * _instFlux / (_instFluxErr * _instFluxErr)); }

    std::shared_ptr<const FittedStar> getFittedStar() const { return _fittedStar; };

    const CcdImage &getCcdImage() const { return *_ccdImage; };

    void setCcdImage(const CcdImage *ccdImage) { _ccdImage = ccdImage; };

    //! Fits may use that to discard outliers
    bool isValid() const { return _valid; }
    //! Fits may use that to discard outliers
    void setValid(bool v) { _valid = v; }
};

/****** MeasuredStarList */

//! A list of MeasuredStar. They are usually filled in Associations::AddImage
class MeasuredStarList : public StarList<MeasuredStar> {
public:
    MeasuredStarList(){};

    void setCcdImage(const CcdImage *_ccdImage);
};

typedef MeasuredStarList::const_iterator MeasuredStarCIterator;
typedef MeasuredStarList::iterator MeasuredStarIterator;

BaseStarList &Measured2Base(MeasuredStarList &This);
BaseStarList *Measured2Base(MeasuredStarList *This);
const BaseStarList &Measured2Base(const MeasuredStarList &This);
const BaseStarList *Measured2Base(const MeasuredStarList *This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_MEASURED_STAR_H
