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

#ifndef LSST_JOINTCAL_PHOTOMETRY_TRANSFORM_H
#define LSST_JOINTCAL_PHOTOMETRY_TRANSFORM_H

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

/**
 * A photometric transform, defined in terms of the input flux or magnitude.
 *
 * Unit agnostic: a higher level Model must keep track of the units going into and out of of its Transforms.
 *
 * @see lsst::afw::image::PhotoCalib
 */
class PhotometryTransform {
public:
    /// Return the transform of value at (x,y).
    virtual double transform(double x, double y, double value) const = 0;

    /// Return the transformed value at Point(x,y).
    double transform(Point const &in, double value) const { return transform(in.x, in.y, value); }

    /// Return the transformed valueErr at Point(x,y).
    virtual double transformError(double x, double y, double value, double valueErr) const = 0;

    /// Return the transformed valueErr at Point(x,y).
    double transformError(Point const &in, double value, double valueErr) const {
        return transformError(in.x, in.y, value, valueErr);
    }

    /// dumps the transform coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    friend std::ostream &operator<<(std::ostream &s, PhotometryTransform const &transform) {
        transform.dump(s);
        return s;
    }

    /// Return the number of parameters (used to compute chisq)
    virtual std::size_t getNpar() const = 0;

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
    virtual std::shared_ptr<PhotometryTransform> clone() const = 0;

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  x        The x coordinate to compute at (in the appropriate units for this transform).
     * @param[in]  y        The y coordinate to compute at (in the appropriate units for this transform).
     * @param[in]  value    The instrument flux or magnitude to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void computeParameterDerivatives(double x, double y, double value,
                                             Eigen::Ref<Eigen::VectorXd> derivatives) const = 0;

    /// Get a copy of the parameters of this model, in the same order as `offsetParams`.
    virtual Eigen::VectorXd getParameters() const = 0;
};

/**
 * Photometry offset independent of position. Abstract class.
 */
class PhotometryTransformSpatiallyInvariant : public PhotometryTransform {
public:
    explicit PhotometryTransformSpatiallyInvariant(double value) : _value(value) {}

    /// @copydoc PhotometryTransform::dump
    void dump(std::ostream &stream = std::cout) const override { stream << std::setprecision(10) << _value; }

    /// @copydoc PhotometryTransform::getNpar
    std::size_t getNpar() const override { return 1; }

    /// @copydoc PhotometryTransform::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override { _value -= delta[0]; };

    /// @copydoc PhotometryTransform::getParameters
    Eigen::VectorXd getParameters() const override {
        Eigen::VectorXd parameters(1);
        parameters[0] = _value;
        return parameters;
    }

protected:
    double getValue() const { return _value; }

private:
    /// value of this transform at all locations.
    double _value;
};

/**
 * Photometric offset independent of position, defined as (fluxMag0)^-1.
 *
 * initialCalibFlux * SpatiallyInvariantTransform -> correctedFlux
 *
 */
class FluxTransformSpatiallyInvariant : public PhotometryTransformSpatiallyInvariant {
public:
    explicit FluxTransformSpatiallyInvariant(double value = 1)
            : PhotometryTransformSpatiallyInvariant(value) {}

    /// @copydoc PhotometryTransform::transform
    double transform(double x, double y, double value) const override { return value * getValue(); }

    /// @copydoc PhotometryTransform::transformError
    double transformError(double x, double y, double value, double valueErr) const override {
        return getValue() * valueErr;
    }

    /// @copydoc PhotometryTransform::clone
    std::shared_ptr<PhotometryTransform> clone() const override {
        return std::make_shared<FluxTransformSpatiallyInvariant>(getValue());
    }

    /// @copydoc PhotometryTransform::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        // the derivative of a spatially constant transform w.r.t. that value is just the value.
        derivatives[0] = value;
    }
};

/**
 * Photometric offset independent of position, defined as -2.5 * log(flux / fluxMag0).
 *
 * initialMagnitude + SpatiallyInvariantTransform -> correctedMagnitude
 *
 */
class MagnitudeTransformSpatiallyInvariant : public PhotometryTransformSpatiallyInvariant {
public:
    explicit MagnitudeTransformSpatiallyInvariant(double value = 0)
            : PhotometryTransformSpatiallyInvariant(value) {}

    /// @copydoc PhotometryTransformSpatiallyInvariant::transform
    double transform(double x, double y, double mag) const override { return mag + getValue(); }

    /// @copydoc PhotometryTransform::transformError
    double transformError(double x, double y, double value, double valueErr) const override {
        return valueErr;
    }

    /// @copydoc PhotometryTransformSpatiallyInvariant::clone
    std::shared_ptr<PhotometryTransform> clone() const override {
        return std::make_shared<MagnitudeTransformSpatiallyInvariant>(getValue());
    }

    /// @copydoc PhotometryTransformSpatiallyInvariant::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        // the derivative of a spatially constant transform w.r.t. that value is 1.
        derivatives[0] = 1;
    }
};

/**
 * nth-order 2d Chebyshev photometry transform.
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
class PhotometryTransformChebyshev : public PhotometryTransform {
public:
    /**
     * Create a Chebyshev transform with terms up to order in (x*y).
     *
     * @param[in]  order  The maximum order in (x*y).
     * @param[in]  bbox   The bounding box it is valid within, to rescale it to [-1,1].
     * @param[in]  identity If true, set a_0,0==1, otherwise all coefficients are 0.
     */
    PhotometryTransformChebyshev(size_t order, geom::Box2D const &bbox, bool identity);

    /**
     * Create a Chebyshev transform with the specified coefficients.
     *
     * The polynomial order is determined from the number of coefficients, taking only the
     * anti-diagonal upper triangle portion of the passed-in coefficients
     *
     * @param      coefficients  The polynomial coefficients.
     * @param[in]  bbox          The bounding box it is valid within, to rescale it to [-1,1].
     */
    PhotometryTransformChebyshev(ndarray::Array<double, 2, 2> const &coefficients, geom::Box2D const &bbox);

    /// @copydoc PhotometryTransform::transformError
    double transformError(double x, double y, double value, double valueErr) const override { return 0; }

    /// @copydoc PhotometryTransform::dump
    void dump(std::ostream &stream = std::cout) const override { stream << _coefficients; }

    /// @copydoc PhotometryTransform::getNpar
    std::size_t getNpar() const override { return _nParameters; }

    /// @copydoc PhotometryTransform::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// Get a copy of the coefficients of the polynomials, as a 2d array (NOTE: layout is [y][x])
    ndarray::Array<double, 2, 2> getCoefficients() const { return ndarray::copy(_coefficients); }

    /// @copydoc PhotometryTransform::getParameters
    Eigen::VectorXd getParameters() const override;

    ndarray::Size getOrder() const { return _order; }

    geom::Box2D getBBox() const { return _bbox; }

    /// Compute the mean of this tranform on the bbox (default to our bbox).
    double mean(geom::Box2D const &bbox) const;
    /// @overload mean(geom::Box2D const &bbox) const;
    double mean() const;

    // Compute the integral of this function over a bounding-box (default to our bbox).
    double integrate(geom::Box2D const &bbox) const;
    /// @overload integrate(geom::Box2D const &bbox) const;
    double integrate() const;

protected:
    /**
     * Return the value of this polynomial at x,y. For use in the sublcass transform() methods.
     */
    double computeChebyshev(double x, double y) const;

    /**
     * Set the derivatives of this polynomial at x,y. For use in the sublcass computeParameterDerivatives()
     * methods.
     */
    void computeChebyshevDerivatives(double x, double y, Eigen::Ref<Eigen::VectorXd> derivatives) const;

private:
    geom::Box2D _bbox;                        // the domain of this function
    geom::AffineTransform _toChebyshevRange;  // maps points from the bbox to [-1,1]x[-1,1]

    ndarray::Array<double, 2, 2> _coefficients;  // shape=(order+1, order+1)
    ndarray::Size _order;
    ndarray::Size _nParameters;

    /// Evaluate one limit of a definite 2-d integral (sum four of these to get the full integral).
    double oneIntegral(double x, double y) const;
};

/**
 * nth-order 2d Chebyshev photometry transform, times the input flux.
 */
class FluxTransformChebyshev : public PhotometryTransformChebyshev {
public:
    FluxTransformChebyshev(size_t order, geom::Box2D const &bbox)
            : PhotometryTransformChebyshev(order, bbox, true) {}

    FluxTransformChebyshev(ndarray::Array<double, 2, 2> const &coefficients, geom::Box2D const &bbox)
            : PhotometryTransformChebyshev(coefficients, bbox) {}

    /// @copydoc PhotometryTransform::transform
    double transform(double x, double y, double value) const override {
        return value * computeChebyshev(x, y);
    }

    /// @copydoc PhotometryTransform::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        computeChebyshevDerivatives(x, y, derivatives);
        derivatives *= value;
    }

    /// @copydoc PhotometryTransform::clone
    std::shared_ptr<PhotometryTransform> clone() const override {
        return std::make_shared<FluxTransformChebyshev>(getCoefficients(), getBBox());
    }
};

/**
 * nth-order 2d Chebyshev photometry transform, plus the input flux.
 */
class MagnitudeTransformChebyshev : public PhotometryTransformChebyshev {
public:
    MagnitudeTransformChebyshev(size_t order, geom::Box2D const &bbox)
            : PhotometryTransformChebyshev(order, bbox, false) {}

    MagnitudeTransformChebyshev(ndarray::Array<double, 2, 2> const &coefficients, geom::Box2D const &bbox)
            : PhotometryTransformChebyshev(coefficients, bbox) {}

    /// @copydoc PhotometryTransform::transform
    double transform(double x, double y, double value) const override {
        return value + computeChebyshev(x, y);
    }

    /// @copydoc PhotometryTransform::computeParameterDerivatives
    void computeParameterDerivatives(double x, double y, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        // The derivatives here are independent of value
        computeChebyshevDerivatives(x, y, derivatives);
    }

    /// @copydoc PhotometryTransform::clone
    std::shared_ptr<PhotometryTransform> clone() const override {
        return std::make_shared<FluxTransformChebyshev>(getCoefficients(), getBBox());
    }
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_TRANSFORM_H
