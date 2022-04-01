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

#ifndef LSST_JOINTCAL_BASE_STAR_H
#define LSST_JOINTCAL_BASE_STAR_H

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <sstream>

#include "lsst/utils/Magnitude.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

namespace {

/// Compute magnitude error from instFlux and instFlux error.
double magErrFromFluxErr(double instFlux, double instFluxErr) {
    return 2.5 / std::log(10) * (instFluxErr / instFlux);
}

}  // namespace

//! The base class for handling stars. Used by all matching routines.
class BaseStar : public FatPoint {
public:
    BaseStar() {
        x = 0;
        y = 0;
        _flux = 0;
    };
    //! constructor
    BaseStar(double xx, double yy, double flux, double fluxErr)
            : FatPoint(xx, yy),
              _flux(flux),
              _fluxErr(fluxErr),
              _mag(utils::nanojanskyToABMagnitude(flux)),
              _magErr(magErrFromFluxErr(flux, fluxErr)){};
    BaseStar(Point const &point, double flux, double fluxErr, double mag, double magErr)
            : FatPoint(point),
              _flux(flux),
              _fluxErr(fluxErr),
              _mag(utils::nanojanskyToABMagnitude(flux)),
              _magErr(magErrFromFluxErr(flux, fluxErr)){};

    //! access stuff.
    double getX() const { return x; }
    //!
    double getY() const { return y; }

    //! allows std::cout << aBaseStar;
    friend std::ostream &operator<<(std::ostream &stream, BaseStar const &s) {
        s.print(stream);
        return stream;
    }

    virtual void print(std::ostream &out) const {
        FatPoint::print(out);
        out << " flux: " << std::setprecision(9) << _flux << " fluxErr: " << _fluxErr;
    }

    BaseStar &operator=(Point const &point) {
        x = point.x;
        y = point.y;
        return (*this);
    };

    static const char *typeName() { return "BaseStar"; }

    virtual ~BaseStar() = default;;

    double getFlux() const { return _flux; }
    double &getFlux() { return _flux; }
    void setFlux(double flux) { _flux = flux; }

    double getFluxErr() const { return _fluxErr; }
    void setFluxErr(double fluxErr) { _fluxErr = fluxErr; }

    double getMag() const { return _mag; }
    double &getMag() { return _mag; }

    double getMagErr() const { return _magErr; }
    void setMagErr(double magErr) { _magErr = magErr; }

protected:
    // on-sky flux, in nanojansky
    double _flux;
    double _fluxErr{};

    // on-sky magnitude
    double _mag{};
    double _magErr{};
};

using BaseStarList = StarList<BaseStar>;

using BaseStarCIterator = BaseStarList::const_iterator;
using BaseStarIterator = BaseStarList::iterator;
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_BASE_STAR_H
