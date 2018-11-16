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

#ifndef LSST_JOINTCAL_POINT_H
#define LSST_JOINTCAL_POINT_H
#include <iostream>
#include <cmath>

namespace lsst {
namespace jointcal {

/*! \file */

//! A point in a plane.
class Point {
public:
    virtual ~Point() {}

    //! coordinate
    double x, y;

    //! - contructor
    Point() : x(0), y(0){};

    //! - contructor
    Point(double xx, double yy) : x(xx), y(yy){};

    //! -
    double Distance(const Point& other) const {
        return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    };

    //! distance squared to other
    double computeDist2(const Point& other) const {
        return ((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    };

    //! Sum
    Point operator+(const Point& Right) const { return Point(x + Right.x, y + Right.y); }

    //! Difference
    Point operator-(const Point& Right) const { return Point(x - Right.x, y - Right.y); }

    //! utility
    virtual void dump(std::ostream& s = std::cout) const { s << "x: " << x << " y: " << y; }

    //! -
    friend std::ostream& operator<<(std::ostream& stream, const Point& point) {
        point.dump(stream);
        return stream;
    }
};

std::ostream& operator<<(std::ostream& stream, const Point& point);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_POINT_H
