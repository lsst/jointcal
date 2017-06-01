// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
#define LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H

#include <iostream>
#include <sstream>

// #include "lsst/jointcal/Point"

namespace lsst {
namespace jointcal {

class PhotometryTransfo {
public:
    /// Apply the transform to (x,y), and puts the result in out
    virtual void apply(double x, double y, double &out) const = 0;

    // double apply(const Point in) const {
    //     double result;
    //     apply(in.x, in.y, result);
    //     return result
    // }

    /// dumps the transfo coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    std::string __str__() { std::stringstream s; dump(s); return s.str(); }

    /// Return the number of parameters (to compute chi2's)
    virtual int getNpar() const {return 0;}
};

/*
 * Constant photmetric offset across an entire ccd.
 */
class ConstantPhotometryTransfo : public PhotometryTransfo {
public:
    int getNpar() const {return 1;}

    void apply(double x, double y, double &out) const {out = value;}

private:
    /// value of this transform at all locations.
    double value;
};

}} // namespaces

#endif // LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
