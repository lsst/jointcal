// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FAT_POINT_H
#define LSST_JOINTCAL_FAT_POINT_H

#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

//! A Point with uncertainties
class FatPoint : public Point {
public:
    double vx, vy, vxy;

    FatPoint() : Point() {
        vx = vy = 1;
        vxy = 0;
    };

    FatPoint(Point const& P, double Vx = 1, double Vy = 1, double Vxy = 0)
            : Point(P), vx(Vx), vy(Vy), vxy(Vxy){};

    FatPoint(const double X, const double Y, const double Vx = 1, const double Vy = 1, const double Vxy = 0)
            : Point(X, Y), vx(Vx), vy(Vy), vxy(Vxy){};

    void dump(std::ostream& s = std::cout) const {
        Point::dump(s);
        s << " vxx,vyy,vxy " << vx << ' ' << vy << ' ' << vxy;
    }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FAT_POINT_H
