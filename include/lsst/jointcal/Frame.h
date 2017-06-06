// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FRAME_H
#define LSST_JOINTCAL_FRAME_H

#include <iostream>
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

typedef enum { WholeSizeFrame, ClippedSizeFrame } WhichFrame;

typedef enum { LargeFrame, SmallFrame } WhichTransformed;

//! rectangle with sides parallel to axes.
/*! when Frame's are used to define subparts of images, xMin and xMax refer
  to the first and last pixels in the subimage */

class Frame {
public:
    //! coordinate of boundary.
    double xMin, xMax, yMin, yMax;

    //! Default constructor
    Frame();

    //! this one is dangerous: you may swap the 2 middle arguments.
    //! Prefer next one
    Frame(double xMin, double yMin, double xMax, double yMax);

    //! typical use: Frame(Point(xmin,ymin),Point(xmax,ymax))
    Frame(const Point &lowerLeft, const Point &upperRight);

    //! size along x axis
    double getWidth() const { return xMax - xMin; }

    //! size along y axis
    double getHeight() const { return yMax - yMin; }

    //! Center  of the frame
    Point getCenter() const { return Point((xMax + xMin) * 0.5, (yMax + yMin) * 0.5); }

    //! intersection of Frame's.
    Frame operator*(const Frame &right) const; /* intersection : a = b n c */

    //! intersection of Frame's
    Frame &operator*=(const Frame &right); /* intersection : a = a n b */

    //! union of Frames
    Frame operator+(const Frame &right) const; /* union : a = b u c */

    //! union of Frames
    Frame &operator+=(const Frame &right); /* intersection : a = a u b */

    //! shrinks the frame (if marginSize>0), enlarges it (if marginSize<0).
    void cutMargin(const double marginSize);

    //! shrinks the frame (if marginSize>0), enlarges it (if marginSize<0).
    void cutMargin(const double marginX, const double marginY);

    //! necessary for comparisons (!= is defined from this one implicitely)
    bool operator==(const Frame &right) const;

    //! comparison
    bool operator!=(const Frame &right) const { return !(*this == right); }

    //! rescale it. The center does not move.
    Frame rescale(const double factor) const;

    // the area.
    double getArea() const;

    //! inside?
    bool inFrame(double x, double y) const;

    //! same as above
    bool inFrame(const Point &point) const { return inFrame(point.x, point.y); }

    //! distance to closest boundary.
    double minDistToEdges(const Point &point) const;

    void dump(std::ostream &stream = std::cout) const;

    //! allows \verbatim std::cout << frame; \endverbatim.
    friend std::ostream &operator<<(std::ostream &stream, const Frame &right) {
        right.dump(stream);
        return stream;
    };

private:
    void order();
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_FRAME_H
