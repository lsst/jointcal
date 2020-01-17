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

    FatPoint(const Point& P, double Vx = 1, double Vy = 1, double Vxy = 0)
            : Point(P), vx(Vx), vy(Vy), vxy(Vxy){};

    FatPoint(const double X, const double Y, const double Vx = 1, const double Vy = 1, const double Vxy = 0)
            : Point(X, Y), vx(Vx), vy(Vy), vxy(Vxy){};

    void print(std::ostream& s = std::cout) const {
        Point::print(s);
        s << " vxx,vyy,vxy " << vx << ' ' << vy << ' ' << vxy;
    }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FAT_POINT_H
