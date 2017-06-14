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
protected:
    // on-sky flux, in Maggies
    double _flux;
    double _fluxErr;

public:
    BaseStar() {
        x = 0;
        y = 0;
        _flux = 0;
    };
    //! constructor
    BaseStar(double xx, double yy, double flux) : FatPoint(xx, yy), _flux(flux){};
    BaseStar(const Point &point, double flux) : FatPoint(point), _flux(flux){};

    //! access stuff.
    double getX() const { return x; }
    //!
    double getY() const { return y; }

    //! allows std::cout << aBaseStar;
    friend std::ostream &operator<<(std::ostream &stream, const BaseStar &s) {
        s.dump(stream);
        return stream;
    }

    virtual std::string __str__() const {
        std::stringstream s;
        dump(s);
        return s.str();
    }

    virtual void dump(std::ostream &stream = std::cout) const {
        stream << "x: " << x << " y: " << y << " flux: " << _flux;
    }

    BaseStar &operator=(const Point &point) {
        this->x = point.x;
        this->y = point.y;
        return (*this);
    };

    static const char *typeName() { return "BaseStar"; }

    virtual ~BaseStar(){};

    double getFlux() const { return _flux; }
    double &getFlux() { return _flux; }
    void setFlux(double flux) { _flux = flux; }

    double getFluxErr() const { return _fluxErr; }
    // double &getFluxErr() { return _fluxErr; }
    void setFluxErr(double fluxErr) { _fluxErr = fluxErr; }
};

//! enables to sort easily a starList (of anything that derives from BaseStar)
bool decreasingFlux(const BaseStar *star1, const BaseStar *star2);

int decodeFormat(const char *formatLine, const char *starName);

typedef StarList<BaseStar> BaseStarList;

typedef BaseStarList::const_iterator BaseStarCIterator;
typedef BaseStarList::iterator BaseStarIterator;
}
}  // namespace lsst

#endif  // LSST_JOINTCAL_BASE_STAR_H
