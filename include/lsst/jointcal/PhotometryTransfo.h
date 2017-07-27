// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
#define LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H

#include <iostream>
#include <sstream>
#include <memory>

#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class Point;

class PhotometryTransfoSpatiallyInvariant;

/*
 * A photometric transform, defined as a scale factor of the input calibration.
 *
 *     initialCalibFlux (Maggies) * transfo(x,y) -> correctedFlux (Maggies)
 *
 * @todo Eventually will be defined the same as PhotoCalib:
 *     instFlux (ADU) -> flux (maggies)
 *
 * @seealso lsst::afw::image::PhotoCalib
 */
class PhotometryTransfo {
public:
    /// Apply the transform to instFlux at (x,y), put result in flux
    virtual void apply(double x, double y, double instFlux, double &flux) const = 0;

    double apply(const Point in, double instFlux) const {
        double flux;
        apply(in.x, in.y, instFlux, flux);
        return flux;
    }

    /// dumps the transfo coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /// Return a string describing this transfo. For the pybind11/python layer.
    std::string __str__() const {
        std::stringstream s;
        dump(s);
        return s.str();
    }

    /// Return the number of parameters (used to compute chisq)
    virtual int getNpar() const { return 0; }

    /// Offset the parameters by some amount during fitting.
    virtual void offsetParams(const double *delta) = 0;

    /// return a copy (allocated by new) of the transformation.
    virtual std::unique_ptr<PhotometryTransfo> clone() const = 0;

    void computeDerivative(const Point &where, PhotometryTransfoSpatiallyInvariant &derivative,
                           const double step = 0.01) const;
};

/*
 * Photometric offset independent of position.
 *
 * initialCalibFlux (Maggies) * SpatiallyInvariantTransfo -> correctedFlux (Maggies)
 *
 * @todo Eventually to be defined as:
 *     instFlux / value = flux
 */
class PhotometryTransfoSpatiallyInvariant : public PhotometryTransfo {
public:
    PhotometryTransfoSpatiallyInvariant(double value = 1) : _value(value) {}

    void apply(double x, double y, double instFlux, double &out) const { out = instFlux * _value; }

    void dump(std::ostream &stream = std::cout) const { stream << _value; }

    int getNpar() const { return 1; }

    void offsetParams(const double *delta) { _value -= *delta; };

    std::unique_ptr<PhotometryTransfo> clone() const {
        return std::unique_ptr<PhotometryTransfo>(new PhotometryTransfoSpatiallyInvariant(_value));
    }

    /// The spatial derivative of a constant zeropoint is 1.
    void computeDerivative(const Point &where, PhotometryTransfoSpatiallyInvariant &derivative,
                           const double step = 0.01) const {
        derivative.setValue(1);
    }

protected:
    void setValue(double value) { _value = value; }

    friend class PhotometryTransfo;

private:
    /// value of this transform at all locations.
    double _value;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
