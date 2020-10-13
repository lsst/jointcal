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

#ifndef LSST_JOINTCAL_MEASURED_STAR_H
#define LSST_JOINTCAL_MEASURED_STAR_H

#include <iostream>

#include "lsst/afw/table/misc.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

namespace {  /// Compute "instrumental magnitude" from instFlux (in counts).
double instMagFromInstFlux(double instFlux) { return -2.5 * std::log10(instFlux); }
}  // namespace

class CcdImage;

/**
 * Sources measured on images.
 *
 * x/y positions and uncertainties (from BaseStar) are in pixels in the image frame.
 * instFlux[Err] are in counts.
 * flux[Err] (from BaseStar) are nJy.
 */
class MeasuredStar : public BaseStar {
public:
    MeasuredStar()
            : BaseStar(),
              _id(0),
              _instFlux(0.),
              _instFluxErr(0.),
              _ccdImage(0),
              _valid(true),
              _xFocal(0.0),
              _yFocal(0.0),
              _instMag(0.),
              _instMagErr(0.) {}

    MeasuredStar(BaseStar const &baseStar)
            : BaseStar(baseStar),
              _id(0),
              _instFlux(0.),
              _instFluxErr(0.),
              _ccdImage(0),
              _valid(true),
              _xFocal(0.0),
              _yFocal(0.0),
              _instMag(0.),
              _instMagErr(0.) {}

    /// No move, allow copy constructor: we may copy the fitted StarLists when associating and
    /// matching catalogs, otherwise Stars should be managed by shared_ptr only.
    MeasuredStar(MeasuredStar const &) = default;
    MeasuredStar(MeasuredStar &&) = delete;
    MeasuredStar &operator=(MeasuredStar const &) = delete;
    MeasuredStar &operator=(MeasuredStar &&) = delete;

    void setFittedStar(std::shared_ptr<FittedStar> fittedStar) {
        if (fittedStar) fittedStar->getMeasurementCount()++;
        _fittedStar = std::move(fittedStar);
    }

    void print(std::ostream &out) const {
        BaseStar::print(out);
        out << " instFlux: " << _instFlux << " instFluxErr: " << _instFluxErr << " id: " << _id
            << " valid: " << _valid;
    }

    void setInstFluxAndErr(double instFlux, double instFluxErr) {
        _instFlux = instFlux;
        _instMag = instMagFromInstFlux(instFlux);
        _instFluxErr = instFluxErr;
        _instMagErr = magErrFromFluxErr(instFlux, instFluxErr);
    }

    double getInstFlux() const { return _instFlux; }
    double getInstFluxErr() const { return _instFluxErr; }
    double getInstMag() const { return _instMag; }
    double getInstMagErr() const { return _instMagErr; }

    void setId(afw::table::RecordId id) { _id = id; }
    afw::table::RecordId getId() { return _id; }

    //! the inverse of the mag variance
    double getMagWeight() const { return (_instFlux * _instFlux / (_instFluxErr * _instFluxErr)); }

    double getXFocal() const { return _xFocal; }
    void setXFocal(double xFocal) { _xFocal = xFocal; }
    double getYFocal() const { return _yFocal; }
    void setYFocal(double yFocal) { _yFocal = yFocal; }

    std::shared_ptr<FittedStar> getFittedStar() const { return _fittedStar; };

    CcdImage const &getCcdImage() const { return *_ccdImage; };

    void setCcdImage(const CcdImage *ccdImage) { _ccdImage = ccdImage; };

    //! Fits may use that to discard outliers
    bool isValid() const { return _valid; }
    //! Fits may use that to discard outliers
    void setValid(bool v) { _valid = v; }

private:
    afw::table::RecordId _id;  // id in original catalog

    // on-chip flux, in ADU
    double _instFlux;
    double _instFluxErr;

    const CcdImage *_ccdImage;
    // Note: _fittedStar is not const, but measuredStar won't modify it.
    std::shared_ptr<FittedStar> _fittedStar;
    bool _valid;

    double _xFocal, _yFocal;

    // on-sensor "magnitudes", used when fitting the MagnitudeModel.
    double _instMag;
    double _instMagErr;
};

/****** MeasuredStarList */

//! A list of MeasuredStar. They are usually filled in Associations::createCcdImage
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
