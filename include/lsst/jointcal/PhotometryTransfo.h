// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
#define LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H

#include <iostream>
#include <sstream>
#include <memory>

#include "Eigen/Core"
#include "ndarray.h"

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class Point;

class PhotometryTransfoSpatiallyInvariant;

/*
 * A photometric transform, defined as a scalar multiple of the input flux.
 *
 * Unit agnostic: a higher level Model must keep track of the units going into and out of of its Transfos.
 *
 *     inputFlux * transfo(x,y) -> correctedFlux
 *
 * @seealso lsst::afw::image::PhotoCalib
 */
class PhotometryTransfo {
public:
    /// Apply the transform to instFlux at (x,y), put result in flux
    virtual double transform(double x, double y, double instFlux) const = 0;

    /// Return the transformed instFlux at (x,y).
    double transform(Point const &in, double instFlux) const { return transform(in.x, in.y, instFlux); }

    /// dumps the transfo coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    friend std::ostream &operator<<(std::ostream &s, PhotometryTransfo const &transfo) {
        transfo.dump(s);
        return s;
    }

    /// Return the number of parameters (used to compute chisq)
    virtual int getNpar() const { return 0; }

    /**
     * Offset the parameters by some (negative) amount during fitting.
     *
     * Equivalent to `flatten(parameters) -= delta`
     *
     * Ordering of delta is the same as the ordering of the derivatives returned from
     * `computeParameterDerivatives`.
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    /// return a copy (allocated by new) of the transformation.
    virtual std::shared_ptr<PhotometryTransfo> clone() const = 0;

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  x        The x coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  y        The y coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  instFlux     The instrument flux to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void computeParameterDerivatives(double x, double y, double instFlux,
                                             Eigen::Ref<Eigen::VectorXd> derivatives) const = 0;

    /// Get a copy of the parameters of this model, in the same order as `offsetParams`.
    virtual Eigen::VectorXd getParameters() const = 0;
};

/*
 * Photometric offset independent of position, defined as (fluxMag0)^-1.
 *
 * initialCalibFlux * SpatiallyInvariantTransfo -> correctedFlux
 *
 */
class PhotometryTransfoSpatiallyInvariant : public PhotometryTransfo {
public:
    PhotometryTransfoSpatiallyInvariant(double value = 1) : _value(value) {}

    /// @copydoc PhotometryTransfo::transform
    double transform(double x, double y, double instFlux) const override { return instFlux * _value; }

    /// @copydoc PhotometryTransfo::dump
    void dump(std::ostream &stream = std::cout) const override { stream << _value; }

    /// @copydoc PhotometryTransfo::getNpar
    int getNpar() const override { return 1; }

    /// @copydoc PhotometryTransfo::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override { _value -= delta[0]; };

    /// @copydoc PhotometryTransfo::clone
    std::shared_ptr<PhotometryTransfo> clone() const override {
        return std::make_shared<PhotometryTransfoSpatiallyInvariant>(_value);
    }

    /// @copydoc PhotometryTransfo::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double instFlux,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        // the derivative of a spatially constant transfo w.r.t. that value is just the instFlux.
        derivatives[0] = instFlux;
    }

    /// @copydoc PhotometryTransfo::getParameters
    Eigen::VectorXd getParameters() const override {
        Eigen::VectorXd parameters(1);
        parameters[0] = _value;
        return parameters;
    }

protected:
    void setValue(double value) { _value = value; }

    friend class PhotometryTransfo;

private:
    /// value of this transform at all locations.
    double _value;
};

/**
 * nth-order 2d Chebyshev photometry transfo.
 *
 * The 2-d Chebyshev polynomial used here is defined as:
 *
 *  @f[
 *  f(x,y) = \sum_i \sum_j a_{i,j} T_i(x) T_j(y)
 *  @f]
 *
 * where @f$T_n(x)@f$ is the n-th order Chebyshev polynomial of @f$x@f$ and
 * @f$a_{i,j}@f$ is the corresponding coefficient of the (i,j) polynomial term.
 *
 * Note that the polynomial order=n means that the highest terms will be of the form:
 *   @f[
 *   a_{0,n}*x^n*y^0, a_{n-1,1}*x^(n-1)*y^1, ..., a_{1,n-1}*x^1*y^(n-1), a_{n,0}*x^0*y^n
 *   @f]
 */
class PhotometryTransfoChebyshev : public PhotometryTransfo {
public:
    /**
     * Create an identity (a_0,0==1) Chebyshev transfo with terms up to order in (x*y).
     *
     * @param[in]  order  The maximum order in (x*y).
     * @param[in]  bbox    The bounding box it is valid within, to rescale it to [-1,1].
     */
    PhotometryTransfoChebyshev(size_t order, afw::geom::Box2D const &bbox);

    /**
     * Create a Chebyshev transfo with the specified coefficients.
     *
     * The polynomial order is determined from the number of coefficients.
     *
     * @param      coefficients  The polynomial coefficients.
     * @param[in]  bbox          The bounding box it is valid within, to rescale it to [-1,1].
     */
    PhotometryTransfoChebyshev(ndarray::Array<double, 2, 2> const &coefficients,
                               afw::geom::Box2D const &bbox);

    /// @copydoc PhotometryTransfo::transform
    double transform(double x, double y, double instFlux) const override;

    /// @copydoc PhotometryTransfo::dump
    void dump(std::ostream &stream = std::cout) const override { stream << _coefficients; }

    /// @copydoc PhotometryTransfo::getNpar
    int getNpar() const override { return _nParameters; }

    /// @copydoc PhotometryTransfo::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryTransfo::clone
    std::shared_ptr<PhotometryTransfo> clone() const override {
        return std::make_shared<PhotometryTransfoChebyshev>(ndarray::copy(_coefficients), _bbox);
    }

    /// @copydoc PhotometryTransfo::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double instFlux,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override;

    /// Get a copy of the coefficients of the polynomials, as a 2d array (NOTE: layout is [y][x])
    ndarray::Array<double, 2, 2> getCoefficients() { return ndarray::copy(_coefficients); }

    /// @copydoc PhotometryTransfo::getParameters
    Eigen::VectorXd getParameters() const override;

    ndarray::Size getOrder() const { return _order; }

    afw::geom::Box2D getBBox() const { return _bbox; }

    // Compute the mean of this function over its bounding-box.
    double mean() const;

private:
    afw::geom::Box2D _bbox;                        // the domain of this function
    afw::geom::AffineTransform _toChebyshevRange;  // maps points from the bbox to [-1,1]x[-1,1]

    ndarray::Array<double, 2, 2> _coefficients;  // shape=(order+1, order+1)
    ndarray::Size _order;
    ndarray::Size _nParameters;

    // Compute the integral of this function over its bounding-box.
    double integrate() const;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
