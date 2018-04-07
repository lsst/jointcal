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

    FatPoint(const Point& p, double vx = 1, double vy = 1, double vxy = 0)
            : Point(p), vx(vx), vy(vy), vxy(vxy){};

    FatPoint(const double x, const double y, const double vx = 1, const double vy = 1, const double vxy = 0)
            : Point(x, y), vx(vx), vy(vy), vxy(vxy){};

    void dump(std::ostream& s = std::cout) const {
        Point::dump(s);
        s << " vxx,vyy,vxy " << vx << ' ' << vy << ' ' << vxy;
    }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FAT_POINT_H
