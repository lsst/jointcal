// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_BASE_STAR_H
#define LSST_JOINTCAL_BASE_STAR_H

#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

#define MEMPIX2DISK 1

#define DECALAGE_IJ_XY 0.
#define DECALAGE_XY_IJ 0.

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
            : FatPoint(xx, yy), _flux(flux), _fluxErr(fluxErr){};
    BaseStar(Point const &point, double flux, double fluxErr)
            : FatPoint(point), _flux(flux), _fluxErr(fluxErr){};

    //! access stuff.
    double getX() const { return x; }
    //!
    double getY() const { return y; }

    //! allows std::cout << aBaseStar;
    friend std::ostream &operator<<(std::ostream &stream, BaseStar const &s) {
        s.dump(stream);
        return stream;
    }

    virtual void dump(std::ostream &stream = std::cout) const {
        stream << "x: " << x << " y: " << y << " flux: " << _flux << " fluxErr: " << _fluxErr;
    }

    BaseStar &operator=(Point const &point) {
        x = point.x;
        y = point.y;
        return (*this);
    };

    static const char *typeName() { return "BaseStar"; }

    virtual ~BaseStar(){};

    double getFlux() const { return _flux; }
    double &getFlux() { return _flux; }
    void setFlux(double flux) { _flux = flux; }

    double getFluxErr() const { return _fluxErr; }
    void setFluxErr(double fluxErr) { _fluxErr = fluxErr; }

protected:
    // on-sky flux, in Maggies
    double _flux;
    double _fluxErr;
};

//! enables to sort easily a starList (of anything that derives from BaseStar)
bool decreasingFlux(BaseStar const *star1, BaseStar const *star2);

int decodeFormat(char const *formatLine, char const *starName);

typedef StarList<BaseStar> BaseStarList;

typedef BaseStarList::const_iterator BaseStarCIterator;
typedef BaseStarList::iterator BaseStarIterator;
}
}  // namespace lsst

#endif  // LSST_JOINTCAL_BASE_STAR_H
