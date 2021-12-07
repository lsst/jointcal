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

#include <iostream>

#include "lsst/jointcal/Frame.h"

namespace lsst {
namespace jointcal {

using namespace std;

/****************** Frame class methods ***********************/
Frame::Frame(const Point &lowerLeft, const Point &upperRight) {
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
double Frame::minDistToEdges(const Point &point) const {
    return min(min(point.x - xMin, xMax - point.x) /* minx */,
               min(point.y - yMin, yMax - point.y) /* miny */);
}

Frame Frame::operator*(const Frame &right) const {
    Frame result = *this;
    result *= right;
    return result;
}

Frame &Frame::operator*=(const Frame &right) {
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

Frame Frame::operator+(const Frame &right) const {
    Frame result = *this;
    result += right;
    return result;
}

Frame &Frame::operator+=(const Frame &right) {
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

bool Frame::operator==(const Frame &right) const {
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

void Frame::print(ostream &stream) const {
    stream << "(xMin,yMin), (xMax,yMax): (" << std::setprecision(12) << xMin << ", " << yMin << "), (" << xMax
           << ", " << yMax << ")";
}
}  // namespace jointcal
}  // namespace lsst
