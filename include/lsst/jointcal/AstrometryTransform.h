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

#ifndef LSST_JOINTCAL_ASTROMETRY_TRANSFORM_H
#define LSST_JOINTCAL_ASTROMETRY_TRANSFORM_H

#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/Frame.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

class StarMatchList;
class Frame;
class AstrometryTransformLinear;

//! a virtual (interface) class for geometric transformations.
/*! We implement here One AstrometryTransform interface class, and actual derived
    classes. Composition in the usual (mathematical) sense is provided
    using compose(), and some classes (e.g. AstrometryTransformLinear)
    handle a * operator.  Generic inversion by iteration exists, but
    it is at least 10 times slower than the corresponding "direct
    transformation". If a transform has an analytical inverse, then
    providing inverseTransform is obviously a very good idea. Before
    resorting to inverseTransform, consider using
    StarMatchList::inverseTransform().  AstrometryTransformLinear::inverted() and
    TanPixelToRaDec::inverted() exist.
    The classes also provide derivation and linear approximation.

*/

class AstrometryTransform {
public:
    //!
    virtual void apply(const double xIn, const double yIn, double &xOut, double &yOut) const = 0;

    //! applies the tranfo to in and writes into out. Is indeed virtual.
    void apply(Point const &in, Point &out) const { apply(in.x, in.y, out.x, out.y); }

    //! All these apply(..) shadow the virtual one in derived classes, unless one writes "using
    //! AstrometryTransform::apply".
    Point apply(Point const &in) const {
        double xout, yout;
        apply(in.x, in.y, xout, yout);
        return Point(xout, yout);
    }

    /**
     * Transform a bounding box, taking either the inscribed or circumscribed box.
     *
     * @param[in] inputframe The frame to be transformed.
     * @param[in] inscribed Return the inscribed (true) or circumscribed (false) box.
     *
     * @return The transformed frame.
     */
    Frame apply(Frame const &inputframe, bool inscribed) const;

    //! dumps the transform coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    std::string __str__() {
        std::stringstream s;
        dump(s);
        return s.str();
    }

    //! fits a transform to a std::list of Point pairs (p1,p2, the Point fields in StarMatch).
    /*! After the fit this(p1) yields approximately p2.
      The returned value is the sum of squared residuals.
      If you want to fit a partial transform (e.g. such that
      this(T1(p1)) = T2(p2), use StarMatchList::applyTransform beforehand. */
    virtual double fit(StarMatchList const &starMatchList) = 0;

    void transformStar(FatPoint &in) const { transformPosAndErrors(in, in); }

    //! returns the local jacobian.
    virtual double getJacobian(Point const &point) const { return getJacobian(point.x, point.y); }

    //! returns a copy (allocated by new) of the transformation.
    virtual std::unique_ptr<AstrometryTransform> clone() const = 0;

    /**
     * Return a reduced composition of newTransform = this(right()), or nullptr if it cannot be reduced.
     *
     * "Reduced" in this context means that they are capable of being merged into a single transform,
     * for example, for two polynomials:
     * @f[
     *     f(x) = 1 + x^2, g(x) = -1 + 3x
     * @f]
     * we would have `h = f.composeAndReduce(g) == 2 - 6x + 9x^2`.
     *
     * To be overloaded by derived classes if they can properly reduce the composition.
     *
     * @param  right  The transform to apply first.
     *
     * @returns The new reduced and composed AstrometryTransform, or nullptr if no such reduction is possible.
     */
    virtual std::unique_ptr<AstrometryTransform> composeAndReduce(AstrometryTransform const &right) const;

    //! returns the local jacobian.
    virtual double getJacobian(const double x, const double y) const;

    /**
     * Computes the local Derivative of a transform, w.r.t. position.
     *
     * Step is used for numerical derivation.
     */
    virtual void computeDerivative(Point const &where, AstrometryTransformLinear &derivative,
                                   const double step = 0.01) const;

    //! linear (local) approximation.
    virtual AstrometryTransformLinear linearApproximation(Point const &where, const double step = 0.01) const;

    virtual void transformPosAndErrors(const FatPoint &in, FatPoint &out) const;

    //! transform errors (represented as double[3] in order V(xx),V(yy),Cov(xy))
    virtual void transformErrors(Point const &where, const double *vIn, double *vOut) const;

    //! returns an inverse transform. Numerical if not overloaded.
    /*! precision and region refer to the "input" side of this,
      and hence to the output side of the returned AstrometryTransform. */
    virtual std::unique_ptr<AstrometryTransform> inverseTransform(const double precision,
                                                                  const Frame &region) const;

    //! params should be at least Npar() long
    void getParams(double *params) const;

    //!
    void offsetParams(Eigen::VectorXd const &delta);

    //!
    virtual double paramRef(const int i) const;

    //!
    virtual double &paramRef(const int i);

    //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar
    //! for y.
    virtual void paramDerivatives(Point const &where, double *dx, double *dy) const;

    //! Rough inverse.
    /*! Stored by the numerical inverter to guess starting point
       for the trials. Just here to enable overloading. */
    virtual std::unique_ptr<AstrometryTransform> roughInverse(const Frame &region) const;

    //! returns the number of parameters (to compute chi2's)
    virtual int getNpar() const { return 0; }

    /**
     * Create an equivalent AST mapping for this transformation, including an analytic inverse if possible.
     *
     * @param   domain The domain of the transform, to help find an inverse.
     *
     * @return  An AST Mapping that represents this transformation.
     */
    virtual std::shared_ptr<ast::Mapping> toAstMap(jointcal::Frame const &domain) const {
        throw std::logic_error("toAstMap is not implemented for this class.");
    }

    void write(const std::string &fileName) const;

    virtual void write(std::ostream &stream) const;

    virtual ~AstrometryTransform(){};
};

/// Delegates to transform.dump()
std::ostream &operator<<(std::ostream &stream, AstrometryTransform const &transform);

/**
 * Returns a pointer to a composition of transforms, representing `left(right())`.
 *
 * Deletion of returned value to be done by caller.
 *
 * If `left->composeAndReduce(right)` returns NULL, build a AstrometryTransformComposition and return it.
 * This routine implements "run-time" compositions. When there is a possible "reduction" (e.g. compositions
 * of polynomials), compose detects it and returns a genuine AstrometryTransform.
 *
 * @returns The composed transform.
 */
std::unique_ptr<AstrometryTransform> compose(AstrometryTransform const &left,
                                             AstrometryTransform const &right);

/*=============================================================*/
//! A do-nothing transformation. It anyway has dummy routines to mimick a AstrometryTransform

class AstrometryTransformIdentity : public AstrometryTransform {
public:
    //! constructor.
    AstrometryTransformIdentity() {}

    //! xOut = xIn; yOut = yIn !
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const override {
        xOut = xIn;
        yOut = yIn;
    }  // to speed up

    double fit(StarMatchList const &starMatchList) override {
        throw pexExcept::TypeError(
                "AstrometryTransformIdentity is the identity transformation: it cannot be fit to anything.");
    }

    /// @copydoc AstrometryTransform::composeAndReduce
    std::unique_ptr<AstrometryTransform> composeAndReduce(AstrometryTransform const &right) const override {
        return right.clone();
    }

    void dump(std::ostream &stream = std::cout) const override { stream << "x' = x\ny' = y" << std::endl; }

    int getNpar() const override { return 0; }

    std::unique_ptr<AstrometryTransform> clone() const override {
        return std::unique_ptr<AstrometryTransform>(new AstrometryTransformIdentity);
    }

    void computeDerivative(Point const &where, AstrometryTransformLinear &derivative,
                           const double step = 0.01) const override;

    //! linear approximation.
    virtual AstrometryTransformLinear linearApproximation(Point const &where,
                                                          const double step = 0.01) const override;

    /// @copydoc AstrometryTransform::toAstMap
    std::shared_ptr<ast::Mapping> toAstMap(jointcal::Frame const &domain) const override;

    void write(std::ostream &s) const override;

    void read(std::istream &s);

    //    ClassDef(AstrometryTransformIdentity,1)
};

/**
 * @overload compose(AstrometryTransform const &, AstrometryTransform const &)
 *
 * @note If instead left is Identity, this method does the correct thing via
 *       AstrometryTransformIdentity::composeAndReduce().
 */
std::unique_ptr<AstrometryTransform> compose(AstrometryTransform const &left,
                                             AstrometryTransformIdentity const &right);

//! Shorthand test to tell if a transform is a simple integer shift
bool isIntegerShift(const AstrometryTransform *transform);

/*====================   AstrometryTransformPolynomial  =======================*/

//! Polynomial transformation class.
class AstrometryTransformPolynomial : public AstrometryTransform {
public:
    /**
     * Default transform : identity for all orders (>=1 ).
     *
     * @param order The highest total power (x+y) of monomials of this polynomial.
     */
    AstrometryTransformPolynomial(const unsigned order = 1);

    //! Constructs a "polynomial image" from an existing transform, over a specified domain
    AstrometryTransformPolynomial(const AstrometryTransform *transform, const Frame &frame, unsigned order,
                                  unsigned nPoint = 1000);

    /**
     * Constructs a polynomial approximation to an afw::geom::TransformPoint2ToPoint2.
     *
     * @param[in] transform The transform to be approximated.
     * @param[in] domain The valid domain of the transform.
     * @param[in] order The polynomial order to use when approximating.
     * @param[in] nSteps The number of sample points per axis (nSteps^2 total points).
     */
    AstrometryTransformPolynomial(std::shared_ptr<afw::geom::TransformPoint2ToPoint2> transform,
                                  jointcal::Frame const &domain, unsigned const order,
                                  unsigned const nSteps = 50);

    /// Sets the polynomial order (the highest sum of exponents of the largest monomial).
    void setOrder(const unsigned order);
    /// Returns the polynomial order.
    unsigned getOrder() const { return _order; }

    using AstrometryTransform::apply;  // to unhide AstrometryTransform::apply(Point const &)

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const override;

    //! specialised analytic routine
    void computeDerivative(Point const &where, AstrometryTransformLinear &derivative,
                           const double step = 0.01) const override;

    //! a mix of apply and Derivative
    virtual void transformPosAndErrors(const FatPoint &in, FatPoint &out) const override;

    //! total number of parameters
    int getNpar() const override { return 2 * _nterms; }

    //! print out of coefficients in a readable form.
    void dump(std::ostream &stream = std::cout) const override;

    //! guess what
    double fit(StarMatchList const &starMatchList) override;

    //! Composition (internal stuff in quadruple precision)
    AstrometryTransformPolynomial operator*(AstrometryTransformPolynomial const &right) const;

    //! Addition
    AstrometryTransformPolynomial operator+(AstrometryTransformPolynomial const &right) const;

    //! Subtraction
    AstrometryTransformPolynomial operator-(AstrometryTransformPolynomial const &right) const;

    using AstrometryTransform::composeAndReduce;  // to unhide
                                                  // AstrometryTransform::composeAndReduce(AstrometryTransform
                                                  // const &)
    /// @copydoc AstrometryTransform::composeAndReduce
    std::unique_ptr<AstrometryTransform> composeAndReduce(AstrometryTransformPolynomial const &right) const;

    std::unique_ptr<AstrometryTransform> clone() const override {
        return std::unique_ptr<AstrometryTransform>(new AstrometryTransformPolynomial(*this));
    }

    //! access to coefficients (read only)
    double coeff(const unsigned powX, const unsigned powY, const unsigned whichCoord) const;

    //! write access
    double &coeff(const unsigned powX, const unsigned powY, const unsigned whichCoord);

    //! read access, zero if beyond order
    double coeffOrZero(const unsigned powX, const unsigned powY, const unsigned whichCoord) const;

    double determinant() const;

    //!
    double paramRef(const int i) const override;

    //!
    double &paramRef(const int i) override;

    //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar
    //! for y.
    void paramDerivatives(Point const &where, double *dx, double *dy) const override;

    /// @copydoc AstrometryTransform::toAstMap
    std::shared_ptr<ast::Mapping> toAstMap(jointcal::Frame const &domain) const override;

    void write(std::ostream &s) const override;
    void read(std::istream &s);

private:
    double computeFit(StarMatchList const &starMatchList, AstrometryTransform const &shiftToCenter,
                      const bool useErrors);

    unsigned _order;              // The highest sum of exponents of the largest monomial.
    unsigned _nterms;             // number of parameters per coordinate
    std::vector<double> _coeffs;  // the actual coefficients
                                  // both polynomials in a single vector to speed up allocation and copies

    /* use std::vector rather than double * to avoid
       writing copy constructor and "operator =".
       Vect would work as well but introduces a dependence
       that can be avoided */

    /* This routine take a double * for the vector because the array can
       then be allocated on the execution stack, which speeds thing
       up. However this uses Variable Length Array (VLA) which is not
       part of C++, but gcc implements it. */
    void computeMonomials(double xIn, double yIn, double *monomial) const;

    /**
     * Return the sparse coefficients matrix that ast::PolyMap requires.
     *
     * @see ast::PolyMap for details of the structure of this matrix.
     */
    ndarray::Array<double, 2, 2> toAstPolyMapCoefficients() const;
};

/**
 * Approximate the inverse by a polynomial, to some precision.
 *
 * @param      forward    Transform to be inverted.
 * @param[in]  domain     The domain of forward.
 * @param[in]  precision  Require that \f$chi2 / (nsteps^2) < precision^2\f$.
 * @param[in]  maxOrder  The maximum order allowed of the inverse polynomial.
 * @param[in]  nSteps     The number of sample points per axis (nSteps^2 total points).
 *
 * @return  A polynomial that best approximates forward.
 */
std::shared_ptr<AstrometryTransformPolynomial> inversePolyTransform(AstrometryTransform const &forward,
                                                                    Frame const &domain,
                                                                    double const precision,
                                                                    int const maxOrder = 9,
                                                                    unsigned const nSteps = 50);

AstrometryTransformLinear normalizeCoordinatesTransform(const Frame &frame);

/*=============================================================*/
//! implements the linear transformations (6 real coefficients).
class AstrometryTransformLinear : public AstrometryTransformPolynomial {
public:
    using AstrometryTransformPolynomial::apply;  // to unhide AstrometryTransform::apply(Point const &)

    //! the default constructor constructs the do-nothing transformation.
    AstrometryTransformLinear() : AstrometryTransformPolynomial(1){};

    //! This triggers an exception if P.getOrder() != 1
    explicit AstrometryTransformLinear(AstrometryTransformPolynomial const &transform);

    //!  enables to combine linear tranformations: T1=T2*T3 is legal.
    AstrometryTransformLinear operator*(AstrometryTransformLinear const &right) const;

    //! returns the inverse: T1 = T2.inverted();
    AstrometryTransformLinear inverted() const;

    // useful?    double jacobian(const double x, const double y) const { return determinant();}

    //!
    void computeDerivative(Point const &where, AstrometryTransformLinear &derivative,
                           const double step = 0.01) const;
    //!
    AstrometryTransformLinear linearApproximation(Point const &where, const double step = 0.01) const;

    //  void dump(std::ostream &stream = std::cout) const;

    // double fit(StarMatchList const &starMatchList);

    //! Construct a AstrometryTransformLinear from parameters
    AstrometryTransformLinear(const double ox, const double oy, const double aa11, const double aa12,
                              const double aa21, const double aa22);

    //! Handy converter:
    AstrometryTransformLinear(AstrometryTransformIdentity const &) : AstrometryTransformPolynomial(1){};

    std::unique_ptr<AstrometryTransform> clone() const {
        return std::unique_ptr<AstrometryTransform>(new AstrometryTransformLinear(*this));
    }

    std::unique_ptr<AstrometryTransform> inverseTransform(const double precision, const Frame &region) const;

    double A11() const { return coeff(1, 0, 0); }
    double A12() const { return coeff(0, 1, 0); }
    double A21() const { return coeff(1, 0, 1); }
    double A22() const { return coeff(0, 1, 1); }
    double Dx() const { return coeff(0, 0, 0); }
    double Dy() const { return coeff(0, 0, 1); }

protected:
    double &a11() { return coeff(1, 0, 0); }
    double &a12() { return coeff(0, 1, 0); }
    double &a21() { return coeff(1, 0, 1); }
    double &a22() { return coeff(0, 1, 1); }
    double &dx() { return coeff(0, 0, 0); }
    double &dy() { return coeff(0, 0, 1); }

    friend class AstrometryTransform;
    friend class AstrometryTransformIdentity;    // for AstrometryTransform::Derivative
    friend class AstrometryTransformPolynomial;  // // for AstrometryTransform::Derivative

private:
    void setOrder(const unsigned order);  // to hide AstrometryTransformPolynomial::setOrder
};

/*=============================================================*/

//! just here to provide a specialized constructor, and fit.
class AstrometryTransformLinearShift : public AstrometryTransformLinear {
public:
    using AstrometryTransform::apply;  // to unhide AstrometryTransform::apply(Point const &)
    //! Add ox and oy.
    AstrometryTransformLinearShift(double ox = 0., double oy = 0.)
            : AstrometryTransformLinear(ox, oy, 1., 0., 0., 1.) {}
    AstrometryTransformLinearShift(Point const &point)
            : AstrometryTransformLinear(point.x, point.y, 1., 0., 0., 1.){};
    double fit(StarMatchList const &starMatchList);

    int getNpar() const { return 2; }
};

/*=============================================================*/
//! just here to provide a specialized constructor, and fit.
class AstrometryTransformLinearRot : public AstrometryTransformLinear {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)

    AstrometryTransformLinearRot() : AstrometryTransformLinear(){};
    AstrometryTransformLinearRot(const double angleRad, const Point *center = nullptr,
                                 const double scaleFactor = 1.0);
    double fit(StarMatchList const &starMatchList);

    int getNpar() const { return 4; }
};

/*=============================================================*/

//! just here to provide specialized constructors. AstrometryTransformLinear fit routine.
class AstrometryTransformLinearScale : public AstrometryTransformLinear {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)
    //!
    AstrometryTransformLinearScale(const double scale = 1)
            : AstrometryTransformLinear(0.0, 0.0, scale, 0., 0., scale){};
    //!
    AstrometryTransformLinearScale(const double scaleX, const double scaleY)
            : AstrometryTransformLinear(0.0, 0.0, scaleX, 0., 0., scaleY){};

    int getNpar() const { return 2; }
};

/**
 * A AstrometryTransform that holds a SkyWcs
 *
 * This is intended to hold the initial estimate for the WCS. It need not be TAN-SIP,
 * nor exactly representable as a FITS WCS.
 *
 * AstrometryTransformSkyWcs does not inherit from BaseTanWcs for two reasons:
 * - There is no need.
 * - It is not clear how to implement all of the BaseTanWcs interface, especially corrections.
 */
class AstrometryTransformSkyWcs : public AstrometryTransform {
public:
    AstrometryTransformSkyWcs(std::shared_ptr<afw::geom::SkyWcs> skyWcs);

    using AstrometryTransform::apply;

    // Input is x, y pixels; output is ICRS RA, Dec in degrees
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const override;

    void dump(std::ostream &stream = std::cout) const override;

    /// Not implemented; throws pex::exceptions::LogicError
    double fit(const StarMatchList &starMatchList) override;

    std::unique_ptr<AstrometryTransform> clone() const override;

    std::shared_ptr<afw::geom::SkyWcs> getSkyWcs() const { return _skyWcs; }

private:
    std::shared_ptr<afw::geom::SkyWcs> _skyWcs;
};

/*==================WCS's transform's =====================================*/

class BaseTanWcs : public AstrometryTransform {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)

    BaseTanWcs(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
               const AstrometryTransformPolynomial *corrections = nullptr);

    BaseTanWcs(const BaseTanWcs &original);

    void operator=(const BaseTanWcs &original);

    /// Transform pixels to ICRS RA, Dec in degrees
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! Get the sky origin (CRVAL in FITS WCS terminology) in degrees
    Point getTangentPoint() const;

    //! The Linear part (corresponding to CD's and CRPIX's)
    AstrometryTransformLinear getLinPart() const;

    //! Get a non-owning pointer to the correction transform polynomial
    const AstrometryTransformPolynomial *getCorr() const { return corr.get(); }

    //! Assign the correction polynomial (what it means is left to derived classes)
    void setCorrections(std::unique_ptr<AstrometryTransformPolynomial> corrections);

    //! Get the pixel origin of the WCS (CRPIX in FITS WCS terminology, but zero-based)
    Point getCrPix() const;

    //! Get a transform from pixels to tangent plane (degrees)
    //! This is a linear transform plus the effects of the correction
    virtual AstrometryTransformPolynomial getPixelToTangentPlane() const = 0;

    //! Transform from pixels to tangent plane (degrees)
    virtual void pixToTangentPlane(double xPixel, double yPixel, double &xTangentPlane,
                                   double &yTangentPlane) const = 0;

    ~BaseTanWcs();

protected:
    AstrometryTransformLinear linPixelToTan;  // transform from pixels to tangent plane (degrees)
                                              // a linear approximation centered at the pixel and sky origins
    std::unique_ptr<AstrometryTransformPolynomial> corr;
    double ra0, dec0;   // sky origin (radians)
    double cos0, sin0;  // cos(dec0), sin(dec0)
};

class TanRaDecToPixel;  // the inverse of TanPixelToRaDec.

/**
 * The transformation that handles pixels to sideral transformations (Gnomonic, possibly with polynomial
 * distortions).
 */
class TanPixelToRaDec : public BaseTanWcs {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)
    //! pixToTan describes the transform from pix to tangent plane (degrees). TangentPoint in degrees.
    //! Corrections are applied between Lin and deprojection parts (as in Swarp).
    TanPixelToRaDec(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
                    const AstrometryTransformPolynomial *corrections = nullptr);

    //! the transformation from pixels to tangent plane (degrees)
    AstrometryTransformPolynomial getPixelToTangentPlane() const;

    //! transforms from pixel space to tangent plane (degrees)
    virtual void pixToTangentPlane(double xPixel, double yPixel, double &xTangentPlane,
                                   double &yTangentPlane) const;

    TanPixelToRaDec();

    //! composition with AstrometryTransformLinear
    TanPixelToRaDec operator*(AstrometryTransformLinear const &right) const;

    using AstrometryTransform::composeAndReduce;  // to unhide
                                                  // AstrometryTransform::composeAndReduce(AstrometryTransform
                                                  // const &)
    /// @copydoc AstrometryTransform::composeAndReduce
    std::unique_ptr<AstrometryTransform> composeAndReduce(AstrometryTransformLinear const &right) const;

    //! approximate inverse : it ignores corrections;
    TanRaDecToPixel inverted() const;

    //! Overload the "generic routine" (available for all AstrometryTransform types
    std::unique_ptr<AstrometryTransform> roughInverse(const Frame &region) const;

    //! Inverse transform: returns a TanRaDecToPixel if there are no corrections, or the iterative solver if
    //! there are.
    std::unique_ptr<AstrometryTransform> inverseTransform(const double precision, const Frame &region) const;

    std::unique_ptr<AstrometryTransform> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(StarMatchList const &starMatchList);
};

//! Implements the (forward) SIP distorsion scheme
class TanSipPixelToRaDec : public BaseTanWcs {
public:
    //! pixToTan describes the transform from pix to tangent plane (degrees). TangentPoint in degrees.
    //! Corrections are applied before Lin.
    TanSipPixelToRaDec(AstrometryTransformLinear const &pixToTan, Point const &tangentPoint,
                       const AstrometryTransformPolynomial *corrections = nullptr);

    //! the transformation from pixels to tangent plane (degrees)
    AstrometryTransformPolynomial getPixelToTangentPlane() const;

    //! transforms from pixel space to tangent plane (degrees)
    virtual void pixToTangentPlane(double xPixel, double yPixel, double &xTangentPlane,
                                   double &yTangentPlane) const;

    TanSipPixelToRaDec();

    //! Inverse transform: returns a TanRaDecToPixel if there are no corrections, or the iterative solver if
    //! there are.
    std::unique_ptr<AstrometryTransform> inverseTransform(const double precision, const Frame &region) const;

    std::unique_ptr<AstrometryTransform> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(StarMatchList const &starMatchList);
};

//! This one is the Tangent Plane (called gnomonic) projection (from celestial sphere to tangent plane)
/*! this transform does not implement corrections, since
   they are defined the other way around (from pixels to sky),
   and not invertible analytically. The inversion of tangent
   point WCS (TanPixelToRaDec) is obtained via inverseTransform().
*/

class TanRaDecToPixel : public AstrometryTransform {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)

    //! assume degrees everywhere.
    TanRaDecToPixel(AstrometryTransformLinear const &tan2Pix, Point const &tangentPoint);

    //!
    TanRaDecToPixel();

    //! The Linear part (corresponding to CD's and CRPIX's)
    AstrometryTransformLinear getLinPart() const;

    //! Resets the projection (or tangent) point
    void setTangentPoint(Point const &tangentPoint);

    //! tangent point coordinates (degrees)
    Point getTangentPoint() const;

    //!
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! transform with analytical derivatives
    void transformPosAndErrors(const FatPoint &in, FatPoint &out) const;

    //! exact typed inverse:
    TanPixelToRaDec inverted() const;

    //! Overload the "generic routine" (available for all AstrometryTransform types
    std::unique_ptr<AstrometryTransform> roughInverse(const Frame &region) const;

    //! Inverse transform: returns a TanPixelToRaDec.
    std::unique_ptr<AstrometryTransform> inverseTransform(const double precision, const Frame &region) const;

    void dump(std::ostream &stream) const;

    std::unique_ptr<AstrometryTransform> clone() const;

    double fit(StarMatchList const &starMatchList);

private:
    double ra0, dec0;  // tangent point (radians)
    double cos0, sin0;
    AstrometryTransformLinear linTan2Pix;  // tangent plane (probably degrees) to pixels
};

//! signature of the user-provided routine that actually does the coordinate transform for UserTransform.
typedef void(AstrometryTransformFun)(const double, const double, double &, double &, const void *);

/**
 * A run-time transform that allows users to define a AstrometryTransform with minimal coding
 * (just the apply routine).
 */
class UserTransform : public AstrometryTransform {
public:
    using AstrometryTransform::apply;  // to unhide apply(const Point&)

    //! the transform routine and extra data that it may need.
    UserTransform(AstrometryTransformFun &fun, const void *userData);

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    void dump(std::ostream &stream = std::cout) const;

    double fit(StarMatchList const &starMatchList);

    std::unique_ptr<AstrometryTransform> clone() const;

private:
    AstrometryTransformFun *_userFun;
    const void *_userData;
};

//! The virtual constructor from a file
std::unique_ptr<AstrometryTransform> astrometryTransformRead(const std::string &fileName);
//! The virtual constructor from a file
std::unique_ptr<AstrometryTransform> astrometryTransformRead(std::istream &s);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTROMETRY_TRANSFORM_H
