// -*- C++ -*-
// $Id: frame.cc,v 1.2 2006/12/22 13:35:41 guy Exp $
//
//
//
// Last Modified: $Date: 2006/12/22 13:35:41 $
// By:            $Author: guy $
//
#include <iostream>

#include "lsst/jointcal/Frame.h"

namespace lsst {
namespace jointcal {

using namespace std;

/****************** Frame class methods ***********************/
Frame::Frame(Point const &lowerLeft, Point const &upperRight) {
    *this = Frame(lowerLeft.x, lowerLeft.y, upperRight.x, upperRight.y);
}

Frame::Frame(double xmin, double ymin, double xmax, double ymax) {
    xMin = min(xmin, xmax);
    xMax = max(xmin, xmax);
    yMin = min(ymin, ymax);
    yMax = max(ymin, ymax);
}

Frame::Frame() { xMin = xMax = yMin = yMax = 0; }

/* positive if inside, negative if outside */
double Frame::minDistToEdges(Point const &point) const {
    return min(min(point.x - xMin, xMax - point.x) /* minx */,
               min(point.y - yMin, yMax - point.y) /* miny */);
}

Frame Frame::operator*(Frame const &right) const {
    Frame result = *this;
    result *= right;
    return result;
}

Frame &Frame::operator*=(Frame const &right) {
    Frame rightCopy = right;
    // make sure that coordinates are properly ordered
    order();
    rightCopy.order();
    xMin = max(xMin, rightCopy.xMin);
    xMax = min(xMax, right.xMax);
    yMin = max(yMin, right.yMin);
    yMax = min(yMax, right.yMax);
    // check for an actual overlap. Why was this check added?
    if (xMin > xMax || yMin > yMax) *this = Frame();
    return *this;
}

Frame Frame::operator+(Frame const &right) const {
    Frame result = *this;
    result += right;
    return result;
}

Frame &Frame::operator+=(Frame const &right) {
    Frame rightCopy = right;
    // make sure that coordinates are properly ordered
    order();
    rightCopy.order();
    xMin = min(xMin, rightCopy.xMin);
    xMax = max(xMax, right.xMax);
    yMin = min(yMin, right.yMin);
    yMax = max(yMax, right.yMax);
    return *this;
}

void Frame::order() {
    if (xMin > xMax) swap(xMin, xMax);
    if (yMin > yMax) swap(yMin, yMax);
}

bool Frame::operator==(Frame const &right) const {
    return ((xMin == right.xMin) && (xMax == right.xMax) && (yMin == right.yMin) && (yMax == right.yMax));
}

void Frame::cutMargin(const double marginX, const double marginY) {
    xMin += marginX;
    yMin += marginY;
    xMax -= marginX;
    yMax -= marginY;
}

void Frame::cutMargin(const double marginSize) { cutMargin(marginSize, marginSize); }

Frame Frame::rescale(const double factor) const {
    double hxsize = fabs(factor * 0.5 * (xMax - xMin));
    double xcenter = 0.5 * (xMax + xMin);
    double hysize = fabs(factor * 0.5 * (yMax - yMin));
    double ycenter = 0.5 * (yMax + yMin);
    return Frame(xcenter - hxsize, ycenter - hysize, xcenter + hxsize, ycenter + hysize);
}

double Frame::getArea() const { return fabs((xMax - xMin) * (yMax - yMin)); }

bool Frame::inFrame(double x, double y) const {
    return ((x <= xMax) && (y <= yMax) && (x >= xMin) && (y >= yMin));
}

void Frame::dump(ostream &stream) const {
    stream << "xmin ymin " << xMin << ' ' << yMin << " xmax ymax " << xMax << ' ' << yMax << endl;
}
}  // namespace jointcal
}  // namespace lsst
