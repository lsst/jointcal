// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_GTRANSFO_H
#define LSST_JOINTCAL_GTRANSFO_H

#include <iostream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/jointcal/FatPoint.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

class StarMatchList;
class Frame;
class GtransfoLin;

//! a virtual (interface) class for geometric transformations.
/*! We implement here One Gtransfo interface class, and actual derived
    classes. Composition in the usual (mathematical) sense is provided
    using gtransfoCompose(), and some classes (e.g. GtransfoLin)
    handle a * operator.  Generic inversion by iteration exists, but
    it is at least 10 times slower than the corresponding "direct
    transformation". If a transfo has an analytical inverse, then
    providing inverseTransfo is obviously a very good idea. Before
    resorting to inverseTransfo, consider using
    StarMatchList::inverseTransfo().  GtransfoLin::invert() and
    TanPix2RaDec::invert() exist.
    The classes also provide derivation and linear approximation.

*/

class Gtransfo {
public:
    //!
    virtual void apply(const double xIn, const double yIn, double &xOut, double &yOut) const = 0;

    //! applies the tranfo to in and writes into out. Is indeed virtual.
    void apply(Point const &in, Point &out) const { apply(in.x, in.y, out.x, out.y); }

    //! All these apply(..) shadow the virtual one in derived classes, unless one writes "using
    //! Gtransfo::apply".
    Point apply(Point const &in) const {
        double xout, yout;
        apply(in.x, in.y, xout, yout);
        return Point(xout, yout);
    }

    //! dumps the transfo coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    std::string __str__() {
        std::stringstream s;
        dump(s);
        return s.str();
    }

    //! fits a transfo to a std::list of Point pairs (p1,p2, the Point fields in StarMatch).
    /*! After the fit this(p1) yields approximately p2.
      The returned value is the sum of squared residuals.
      If you want to fit a partial transfo (e.g. such that
      this(T1(p1)) = T2(p2), use StarMatchList::applyTransfo beforehand. */
    virtual double fit(StarMatchList const &starMatchList) = 0;

    //! allows to write MyTransfo(MyStar)
    void transformStar(FatPoint &in) const { transformPosAndErrors(in, in); }

    //! returns the local jacobian.
    virtual double getJacobian(Point const &point) const { return getJacobian(point.x, point.y); }

    //! returns a copy (allocated by new) of the transformation.
    virtual std::unique_ptr<Gtransfo> clone() const = 0;

    //! to be overloaded by derived classes if they can really "reduce" the composition (e.g. composition of
    //! Polynomial can be reduced)
    virtual std::unique_ptr<Gtransfo> reduceCompo(Gtransfo const *right) const;

    //! returns the local jacobian.
    virtual double getJacobian(const double x, const double y) const;

    /**
     * Computes the local Derivative of a transfo, w.r.t. position.
     *
     * Step is used for numerical derivation.
     */
    virtual void computeDerivative(Point const &where, GtransfoLin &derivative,
                                   const double step = 0.01) const;

    //! linear (local) approximation.
    virtual GtransfoLin linearApproximation(Point const &where, const double step = 0.01) const;

    virtual void transformPosAndErrors(FatPoint const &in, FatPoint &out) const;

    //! transform errors (represented as double[3] in order V(xx),V(yy),Cov(xy))
    virtual void transformErrors(Point const &where, double const *vIn, double *vOut) const;

    //! returns an inverse transfo. Numerical if not overloaded.
    /*! precision and region refer to the "input" side of this,
      and hence to the output side of the returned Gtransfo. */
    virtual std::unique_ptr<Gtransfo> inverseTransfo(const double precision, Frame const &region) const;

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
    virtual std::unique_ptr<Gtransfo> roughInverse(Frame const &region) const;

    //! returns the number of parameters (to compute chi2's)
    virtual int getNpar() const { return 0; }

    void write(const std::string &fileName) const;

    virtual void write(std::ostream &stream) const;

    virtual ~Gtransfo(){};
};

//! allows 'stream << Transfo;' (by calling gtransfo.dump(stream)).
std::ostream &operator<<(std::ostream &stream, Gtransfo const &gtransfo);

//! Returns a pointer to a composition. if left->reduceCompo(right) return NULL, builds a GtransfoComposition
//! and returns it. deletion of returned value to be done by caller

std::unique_ptr<Gtransfo> gtransfoCompose(Gtransfo const *left, Gtransfo const *right);

/*=============================================================*/
//! A do-nothing transformation. It anyway has dummy routines to mimick a Gtransfo

class GtransfoIdentity : public Gtransfo {
public:
    //! constructor.
    GtransfoIdentity() {}

    //! xOut = xIn; yOut = yIn !
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const {
        xOut = xIn;
        yOut = yIn;
    };  // to speed up

    double fit(StarMatchList const &starMatchList) {
        throw pexExcept::TypeError(
                "GtransfoIdentity is the identity transformation: it cannot be fit to anything.");
    }

    std::unique_ptr<Gtransfo> reduceCompo(Gtransfo const *right) const { return right->clone(); }
    void dump(std::ostream &stream = std::cout) const { stream << "x' = x\ny' = y" << std::endl; }

    int getNpar() const { return 0; }
    std::unique_ptr<Gtransfo> clone() const { return std::unique_ptr<Gtransfo>(new GtransfoIdentity); }

    void computeDerivative(Point const &where, GtransfoLin &derivative, const double step = 0.01) const;

    //! linear approximation.
    virtual GtransfoLin linearApproximation(Point const &where, const double step = 0.01) const;

    void write(std::ostream &s) const;

    void read(std::istream &s);

    //    ClassDef(GtransfoIdentity,1)
};

//! Shorthand test to tell if a transfo belongs to the GtransfoIdentity class.
bool isIdentity(Gtransfo const *gtransfo);

//! Shorthand test to tell if a transfo is a simple integer shift
bool isIntegerShift(Gtransfo const *gtransfo);

/*====================   GtransfoPoly  =======================*/

//! Polynomial transformation class.
class GtransfoPoly : public Gtransfo {
public:
    //! Default transfo : identity for all degrees (>=1 ). The degree refers to the highest total power (x+y)
    //! of monomials.
    GtransfoPoly(const unsigned degree = 1);

    //! Constructs a "polynomial image" from an existing transfo, over a specified domain
    GtransfoPoly(Gtransfo const *gtransfo, Frame const &frame, unsigned degree, unsigned nPoint = 1000);

    // sets the polynomial degree.
    void setDegree(const unsigned degree);

    using Gtransfo::apply;  // to unhide Gtransfo::apply(Point const &)

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! specialised analytic routine
    void computeDerivative(Point const &where, GtransfoLin &derivative, const double step = 0.01) const;

    //! a mix of apply and Derivative
    virtual void transformPosAndErrors(FatPoint const &in, FatPoint &out) const;

    //! returns degree
    unsigned getDegree() const { return _degree; }

    //! total number of parameters
    int getNpar() const { return 2 * _nterms; }

    //! print out of coefficients in a readable form.
    void dump(std::ostream &stream = std::cout) const;

    //! guess what
    double fit(StarMatchList const &starMatchList);

    //! Composition (internal stuff in quadruple precision)
    GtransfoPoly operator*(GtransfoPoly const &right) const;

    //! Addition
    GtransfoPoly operator+(GtransfoPoly const &right) const;

    //! Subtraction
    GtransfoPoly operator-(GtransfoPoly const &right) const;

    std::unique_ptr<Gtransfo> reduceCompo(Gtransfo const *right) const;

    std::unique_ptr<Gtransfo> clone() const { return std::unique_ptr<Gtransfo>(new GtransfoPoly(*this)); }

    //! access to coefficients (read only)
    double coeff(const unsigned powX, const unsigned powY, const unsigned whichCoord) const;

    //! write access
    double &coeff(const unsigned powX, const unsigned powY, const unsigned whichCoord);

    //! read access, zero if beyond degree
    double coeffOrZero(const unsigned powX, const unsigned powY, const unsigned whichCoord) const;

    double determinant() const;

    //!
    double paramRef(const int i) const;

    //!
    double &paramRef(const int i);

    //! Derivative w.r.t parameters. Derivatives should be al least 2*NPar long. first Npar, for x, last Npar
    //! for y.
    void paramDerivatives(Point const &where, double *dx, double *dy) const;

    void write(std::ostream &s) const;
    void read(std::istream &s);

private:
    double computeFit(StarMatchList const &starMatchList, Gtransfo const &InTransfo, const bool UseErrors);

    unsigned _degree;             // the degree
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
};

//! approximates the inverse by a polynomial, up to required precision.
std::unique_ptr<GtransfoPoly> inversePolyTransfo(Gtransfo const &Direct, Frame const &frame,
                                                 const double Prec);

GtransfoLin normalizeCoordinatesTransfo(Frame const &frame);

/*=============================================================*/
//! implements the linear transformations (6 real coefficients).
class GtransfoLin : public GtransfoPoly {
public:
    using GtransfoPoly::apply;  // to unhide Gtransfo::apply(Point const &)

    //! the default constructor constructs the do-nothing transformation.
    GtransfoLin() : GtransfoPoly(1){};

    //! This triggers an exception if P.degree() != 1
    explicit GtransfoLin(GtransfoPoly const &gtransfoPoly);

    //!  enables to combine linear tranformations: T1=T2*T3 is legal.
    GtransfoLin operator*(GtransfoLin const &right) const;

    //! returns the inverse: T1 = T2.invert();
    GtransfoLin invert() const;

    // useful?    double jacobian(const double x, const double y) const { return determinant();}

    //!
    void computeDerivative(Point const &where, GtransfoLin &derivative, const double step = 0.01) const;
    //!
    GtransfoLin linearApproximation(Point const &where, const double step = 0.01) const;

    //  void dump(std::ostream &stream = std::cout) const;

    // double fit(StarMatchList const &starMatchList);

    //! Construct a GtransfoLin from parameters
    GtransfoLin(const double ox, const double oy, const double aa11, const double aa12, const double aa21,
                const double aa22);

    //! Handy converter:
    GtransfoLin(GtransfoIdentity const &) : GtransfoPoly(1){};

    std::unique_ptr<Gtransfo> clone() const { return std::unique_ptr<Gtransfo>(new GtransfoLin(*this)); }

    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, Frame const &region) const;

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

    friend class Gtransfo;
    friend class GtransfoIdentity;  // for Gtransfo::Derivative
    friend class GtransfoPoly;      // // for Gtransfo::Derivative

private:
    void setDegree(const unsigned degree);  // to hide GtransfoPoly::setDegree
};

/*=============================================================*/

//! just here to provide a specialized constructor, and fit.
class GtransfoLinShift : public GtransfoLin {
public:
    using Gtransfo::apply;  // to unhide Gtransfo::apply(Point const &)
    //! Add ox and oy.
    GtransfoLinShift(double ox = 0., double oy = 0.) : GtransfoLin(ox, oy, 1., 0., 0., 1.) {}
    GtransfoLinShift(Point const &point) : GtransfoLin(point.x, point.y, 1., 0., 0., 1.){};
    double fit(StarMatchList const &starMatchList);

    int getNpar() const { return 2; }
};

/*=============================================================*/
//! just here to provide a specialized constructor, and fit.
class GtransfoLinRot : public GtransfoLin {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)

    GtransfoLinRot() : GtransfoLin(){};
    GtransfoLinRot(const double angleRad, Point const *center = nullptr, const double scaleFactor = 1.0);
    double fit(StarMatchList const &starMatchList);

    int getNpar() const { return 4; }
};

/*=============================================================*/

//! just here to provide specialized constructors. GtransfoLin fit routine.
class GtransfoLinScale : public GtransfoLin {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)
    //!
    GtransfoLinScale(const double scale = 1) : GtransfoLin(0.0, 0.0, scale, 0., 0., scale){};
    //!
    GtransfoLinScale(const double scaleX, const double scaleY)
            : GtransfoLin(0.0, 0.0, scaleX, 0., 0., scaleY){};

    int getNpar() const { return 2; }
};

/**
 * A Gtransfo that holds a SkyWcs
 *
 * This is intended to hold the initial estimate for the WCS. It need not be TAN-SIP,
 * nor exactly representable as a FITS WCS.
 *
 * GtransfoSkyWcs does not inherit from BaseTanWcs for two reasons:
 * - There is no need.
 * - It is not clear how to implement all of the BaseTanWcs interface, especially corrections.
 */
class GtransfoSkyWcs : public Gtransfo {
public:
    GtransfoSkyWcs(std::shared_ptr<afw::geom::SkyWcs> skyWcs);

    using Gtransfo::apply;

    // Input is x, y pixels; output is ICRS RA, Dec in degrees
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const override;

    void dump(std::ostream &stream = std::cout) const override;

    /// Not implemented; throws pex::exceptions::LogicError
    double fit(StarMatchList const &starMatchList) override;

    std::unique_ptr<Gtransfo> clone() const override;

    std::shared_ptr<afw::geom::SkyWcs> getSkyWcs() const { return _skyWcs; }

private:
    std::shared_ptr<afw::geom::SkyWcs> _skyWcs;
};

/*==================WCS's transfo's =====================================*/

class BaseTanWcs : public Gtransfo {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)

    BaseTanWcs(GtransfoLin const &pix2Tan, Point const &tangentPoint,
               GtransfoPoly const *corrections = nullptr);

    BaseTanWcs(BaseTanWcs const &original);

    void operator=(BaseTanWcs const &original);

    /// Transform pixels to ICRS RA, Dec in degrees
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! Get the sky origin (CRVAL in FITS WCS terminology) in degrees
    Point getTangentPoint() const;

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin getLinPart() const;

    //! Get a non-owning pointer to the correction transform polynomial
    GtransfoPoly const *getCorr() const { return corr.get(); }

    //! Assign the correction polynomial (what it means is left to derived classes)
    void setCorrections(std::unique_ptr<GtransfoPoly> corrections);

    //! Get the pixel origin of the WCS (CRPIX in FITS WCS terminology, but zero-based)
    Point getCrPix() const;

    //! Get a transform from pixels to tangent plane (degrees)
    //! This is a linear transform plus the effects of the correction
    virtual GtransfoPoly getPix2TangentPlane() const = 0;

    //! Transform from pixels to tangent plane (degrees)
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const = 0;

    ~BaseTanWcs();

protected:
    GtransfoLin linPix2Tan;  // transform from pixels to tangent plane (degrees)
                             // a linear approximation centered at the pixel and sky origins
    std::unique_ptr<GtransfoPoly> corr;
    double ra0, dec0;   // sky origin (radians)
    double cos0, sin0;  // cos(dec0), sin(dec0)
};

class TanRaDec2Pix;  // the inverse of TanPix2RaDec.

//! the transformation that handles pix to sideral transfos (Gnomonic, possibly with polynomial distortions).
class TanPix2RaDec : public BaseTanWcs {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)
    //! pix2Tan describes the transfo from pix to tangent plane (degrees). TangentPoint in degrees.
    //! Corrections are applied between Lin and deprojection parts (as in Swarp).
    TanPix2RaDec(GtransfoLin const &pix2Tan, Point const &tangentPoint,
                 GtransfoPoly const *corrections = nullptr);

    //! the transformation from pixels to tangent plane (degrees)
    GtransfoPoly getPix2TangentPlane() const;

    //! transforms from pixel space to tangent plane (degrees)
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const;

    TanPix2RaDec();

    //! composition with GtransfoLin
    TanPix2RaDec operator*(GtransfoLin const &right) const;

    std::unique_ptr<Gtransfo> reduceCompo(Gtransfo const *right) const;

    //! approximate inverse : it ignores corrections;
    TanRaDec2Pix invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    std::unique_ptr<Gtransfo> roughInverse(Frame const &region) const;

    //! Inverse transfo: returns a TanRaDec2Pix if there are no corrections, or the iterative solver if there
    //! are.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, Frame const &region) const;

    std::unique_ptr<Gtransfo> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(StarMatchList const &starMatchList);
};

//! Implements the (forward) SIP distorsion scheme
class TanSipPix2RaDec : public BaseTanWcs {
public:
    //! pix2Tan describes the transfo from pix to tangent plane (degrees). TangentPoint in degrees.
    //! Corrections are applied before Lin.
    TanSipPix2RaDec(GtransfoLin const &pix2Tan, Point const &tangentPoint,
                    GtransfoPoly const *corrections = nullptr);

    //! the transformation from pixels to tangent plane (degrees)
    GtransfoPoly getPix2TangentPlane() const;

    //! transforms from pixel space to tangent plane (degrees)
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const;

    TanSipPix2RaDec();

    //! Inverse transfo: returns a TanRaDec2Pix if there are no corrections, or the iterative solver if there
    //! are.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, Frame const &region) const;

    std::unique_ptr<Gtransfo> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(StarMatchList const &starMatchList);
};

//! This one is the Tangent Plane (called gnomonic) projection (from celestial sphere to tangent plane)
/*! this transfo does not implement corrections, since
   they are defined the other way around (from pixels to sky),
   and not invertible analytically. The inversion of tangent
   point WCS (TanPix2RaDec) is obtained via inverseTransfo().
*/

class TanRaDec2Pix : public Gtransfo {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)

    //! assume degrees everywhere.
    TanRaDec2Pix(GtransfoLin const &tan2Pix, Point const &tangentPoint);

    //!
    TanRaDec2Pix();

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin getLinPart() const;

    //! Resets the projection (or tangent) point
    void setTangentPoint(Point const &tangentPoint);

    //! tangent point coordinates (degrees)
    Point getTangentPoint() const;

    //!
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! transform with analytical derivatives
    void transformPosAndErrors(FatPoint const &in, FatPoint &out) const;

    //! exact typed inverse:
    TanPix2RaDec invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    std::unique_ptr<Gtransfo> roughInverse(Frame const &region) const;

    //! Inverse transfo: returns a TanPix2RaDec.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, Frame const &region) const;

    void dump(std::ostream &stream) const;

    std::unique_ptr<Gtransfo> clone() const;

    double fit(StarMatchList const &starMatchList);

private:
    double ra0, dec0;  // tangent point (radians)
    double cos0, sin0;
    GtransfoLin linTan2Pix;  // tangent plane (probably degrees) to pixels
};

//! signature of the user-provided routine that actually does the coordinate transfo for UserTransfo.
typedef void(GtransfoFun)(const double, const double, double &, double &, void const *);

//! a run-time transfo that allows users to define a Gtransfo with minimal coding (just the transfo routine).
class UserTransfo : public Gtransfo {
public:
    using Gtransfo::apply;  // to unhide apply(Point const &)

    //! the transfo routine and extra data that it may need.
    UserTransfo(GtransfoFun &fun, void const *userData);

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    void dump(std::ostream &stream = std::cout) const;

    double fit(StarMatchList const &starMatchList);

    std::unique_ptr<Gtransfo> clone() const;

private:
    GtransfoFun *_userFun;
    void const *_userData;
};

//! The virtual constructor from a file
std::unique_ptr<Gtransfo> gtransfoRead(const std::string &fileName);
//! The virtual constructor from a file
std::unique_ptr<Gtransfo> gtransfoRead(std::istream &s);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_GTRANSFO_H
