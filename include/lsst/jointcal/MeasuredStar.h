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
    unsigned _id;  // id in original catalog

    // on-chip flux, in ADU
    double _instFlux;
    double _instFluxErr;

    const CcdImage *_ccdImage;
    std::shared_ptr<const FittedStar> _fittedStar;
    bool _valid;

    double _xFocal, _yFocal;

public:
    //!
    MeasuredStar()
            : BaseStar(),
              mag(0.),
              wmag(0.),
              _id(0),
              _instFlux(0.),
              _instFluxErr(0.),
              _ccdImage(0),
              _valid(true),
              _xFocal(0.0),
              _yFocal(0.0) {}

    MeasuredStar(BaseStar const &baseStar)
            : BaseStar(baseStar),
              mag(0.),
              wmag(0.),
              _id(0),
              _instFluxErr(0.),
              _ccdImage(0),
              _valid(true),
              _xFocal(0.0),
              _yFocal(0.0) {}

    /// No move, allow copy constructor: we may copy the fitted StarLists when associating and matching
    /// catalogs, otherwise Stars should be managed by shared_ptr only.
    MeasuredStar(MeasuredStar const &) = default;
    MeasuredStar(MeasuredStar &&) = delete;
    MeasuredStar &operator=(MeasuredStar const &) = delete;
    MeasuredStar &operator=(MeasuredStar &&) = delete;

    void setFittedStar(std::shared_ptr<FittedStar> fittedStar) {
        if (fittedStar) fittedStar->getMeasurementCount()++;
        _fittedStar = std::move(fittedStar);
    }

    void dump(std::ostream &stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " instFlux: " << _instFlux << " instFluxErr: " << _instFluxErr << " id: " << _id
               << " valid: " << _valid;
    }

    void setInstFlux(double instFlux) { _instFlux = instFlux; }
    void setInstFluxErr(double instFluxErr) { _instFluxErr = instFluxErr; }

    double getInstFlux() const { return _instFlux; }
    double getInstFluxErr() const { return _instFluxErr; }
    double getMag() const { return mag; }

    void setId(unsigned id) { _id = id; }
    unsigned getId() { return _id; }

    //! the inverse of the mag variance
    double getMagWeight() const { return (_instFlux * _instFlux / (_instFluxErr * _instFluxErr)); }

    double getXFocal() const { return _xFocal; }
    void setXFocal(double xFocal) { _xFocal = xFocal; }
    double getYFocal() const { return _yFocal; }
    void setYFocal(double yFocal) { _yFocal = yFocal; }

    std::shared_ptr<const FittedStar> getFittedStar() const { return _fittedStar; };

    CcdImage const &getCcdImage() const { return *_ccdImage; };

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
