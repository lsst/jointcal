// -*- LSST-C++ -*-
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
    virtual void dump(std::ostream& s = std::cout) const { s << " x " << x << " y " << y; }

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
