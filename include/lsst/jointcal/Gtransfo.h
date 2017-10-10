// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_GTRANSFO_H
#define LSST_JOINTCAL_GTRANSFO_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "Eigen/Core"

#include "lsst/pex/exceptions.h"
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
    void apply(const Point &in, Point &out) const { apply(in.x, in.y, out.x, out.y); }

    //! All these apply(..) shadow the virtual one in derived classes, unless one writes "using
    //! Gtransfo::apply".
    Point apply(const Point &in) const {
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
    virtual double fit(const StarMatchList &starMatchList) = 0;

    //! allows to write MyTransfo(MyStar)
    void transformStar(FatPoint &in) const { transformPosAndErrors(in, in); }

    //! returns the local jacobian.
    virtual double getJacobian(const Point &point) const { return getJacobian(point.x, point.y); }

    //! returns a copy (allocated by new) of the transformation.
    virtual std::unique_ptr<Gtransfo> clone() const = 0;

    //! to be overloaded by derived classes if they can really "reduce" the composition (e.g. composition of
    //! Polynomial can be reduced)
    virtual std::unique_ptr<Gtransfo> reduceCompo(const Gtransfo *right) const;

    //! returns the local jacobian.
    virtual double getJacobian(const double x, const double y) const;

    /**
     * Computes the local Derivative of a transfo, w.r.t. position.
     *
     * Step is used for numerical derivation.
     */
    virtual void computeDerivative(const Point &where, GtransfoLin &derivative,
                                   const double step = 0.01) const;

    //! linear (local) approximation.
    virtual GtransfoLin linearApproximation(const Point &where, const double step = 0.01) const;

    virtual void transformPosAndErrors(const FatPoint &in, FatPoint &out) const;

    //! transform errors (represented as double[3] in order V(xx),V(yy),Cov(xy))
    virtual void transformErrors(const Point &where, const double *vIn, double *vOut) const;

    //! returns an inverse transfo. Numerical if not overloaded.
    /*! precision and region refer to the "input" side of this,
      and hence to the output side of the returned Gtransfo. */
    virtual std::unique_ptr<Gtransfo> inverseTransfo(const double precision, const Frame &region) const;

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
    virtual void paramDerivatives(const Point &where, double *dx, double *dy) const;

    //! Rough inverse.
    /*! Stored by the numerical inverter to guess starting point
       for the trials. Just here to enable overloading. */
    virtual std::unique_ptr<Gtransfo> roughInverse(const Frame &region) const;

    //! returns the number of parameters (to compute chi2's)
    virtual int getNpar() const { return 0; }

    void write(const std::string &fileName) const;

    virtual void write(std::ostream &stream) const;

    virtual ~Gtransfo(){};
};

//! allows 'stream << Transfo;' (by calling gtransfo.dump(stream)).
std::ostream &operator<<(std::ostream &stream, const Gtransfo &gtransfo);

//! Returns a pointer to a composition. if left->reduceCompo(right) return NULL, builds a GtransfoComposition
//! and returns it. deletion of returned value to be done by caller

std::unique_ptr<Gtransfo> gtransfoCompose(const Gtransfo *left, const Gtransfo *right);

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

    double fit(const StarMatchList &starMatchList) {
        throw pexExcept::TypeError(
                "GtransfoIdentity is the identity transformation: it cannot be fit to anything.");
    }

    std::unique_ptr<Gtransfo> reduceCompo(const Gtransfo *right) const { return right->clone(); }
    void dump(std::ostream &stream = std::cout) const { stream << "x' = x\ny' = y" << std::endl; }

    int getNpar() const { return 0; }
    std::unique_ptr<Gtransfo> clone() const { return std::unique_ptr<Gtransfo>(new GtransfoIdentity); }

    void computeDerivative(const Point &where, GtransfoLin &derivative, const double step = 0.01) const;

    //! linear approximation.
    virtual GtransfoLin linearApproximation(const Point &where, const double step = 0.01) const;

    void write(std::ostream &s) const;

    void read(std::istream &s);

    //    ClassDef(GtransfoIdentity,1)
};

//! Shorthand test to tell if a transfo belongs to the GtransfoIdentity class.
bool isIdentity(const Gtransfo *gtransfo);

//! Shorthand test to tell if a transfo is a simple integer shift
bool isIntegerShift(const Gtransfo *gtransfo);

/*====================   GtransfoPoly  =======================*/

//! Polynomial transformation class.
class GtransfoPoly : public Gtransfo {
public:
    using Gtransfo::apply;  // to unhide Gtransfo::apply(const Point &)

private:
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

public:
    //! Default transfo : identity for all degrees (>=1 ). The degree refers to the highest total power (x+y)
    //! of monomials.
    GtransfoPoly(const unsigned degree = 1);

    //! Constructs a "polynomial image" from an existing transfo, over a specified domain
    GtransfoPoly(const Gtransfo *gtransfo, const Frame &frame, unsigned degree, unsigned nPoint = 1000);

    // sets the polynomial degree.
    void setDegree(const unsigned degree);

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! specialised analytic routine
    void computeDerivative(const Point &where, GtransfoLin &derivative, const double step = 0.01) const;

    //! a mix of apply and Derivative
    virtual void transformPosAndErrors(const FatPoint &in, FatPoint &out) const;

    //! returns degree
    unsigned getDegree() const { return _degree; }

    //! total number of parameters
    int getNpar() const { return 2 * _nterms; }

    //! print out of coefficients in a readable form.
    void dump(std::ostream &stream = std::cout) const;

    //! guess what
    double fit(const StarMatchList &starMatchList);

    //! Composition (internal stuff in quadruple precision)
    GtransfoPoly operator*(const GtransfoPoly &right) const;

    //! Addition
    GtransfoPoly operator+(const GtransfoPoly &right) const;

    //! Subtraction
    GtransfoPoly operator-(const GtransfoPoly &right) const;

    std::unique_ptr<Gtransfo> reduceCompo(const Gtransfo *right) const;

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
    void paramDerivatives(const Point &where, double *dx, double *dy) const;

    void write(std::ostream &s) const;
    void read(std::istream &s);

private:
    double computeFit(const StarMatchList &starMatchList, const Gtransfo &InTransfo, const bool UseErrors);
};

//! approximates the inverse by a polynomial, up to required precision.
std::unique_ptr<GtransfoPoly> inversePolyTransfo(const Gtransfo &Direct, const Frame &frame,
                                                 const double Prec);

GtransfoLin normalizeCoordinatesTransfo(const Frame &frame);

/*=============================================================*/
//! implements the linear transformations (6 real coefficients).
class GtransfoLin : public GtransfoPoly {
public:
    using GtransfoPoly::apply;  // to unhide Gtransfo::apply(const Point &)

    //! the default constructor constructs the do-nothing transformation.
    GtransfoLin() : GtransfoPoly(1){};

    //! This triggers an exception if P.degree() != 1
    explicit GtransfoLin(const GtransfoPoly &gtransfoPoly);

    //!  enables to combine linear tranformations: T1=T2*T3 is legal.
    GtransfoLin operator*(const GtransfoLin &right) const;

    //! returns the inverse: T1 = T2.invert();
    GtransfoLin invert() const;

    // useful?    double jacobian(const double x, const double y) const { return determinant();}

    //!
    void computeDerivative(const Point &where, GtransfoLin &derivative, const double step = 0.01) const;
    //!
    GtransfoLin linearApproximation(const Point &where, const double step = 0.01) const;

    //  void dump(std::ostream &stream = std::cout) const;

    // double fit(const StarMatchList &starMatchList);

    //! the constructor that enables to set all parameters independently. Not very useful.
    GtransfoLin(const double ox, const double oy, const double aa11, const double aa12, const double aa21,
                const double aa22);

    //! Handy converter:
    GtransfoLin(const GtransfoIdentity &) : GtransfoPoly(1){};

    std::unique_ptr<Gtransfo> clone() const { return std::unique_ptr<Gtransfo>(new GtransfoLin(*this)); }

    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, const Frame &region) const;

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
    using Gtransfo::apply;  // to unhide Gtransfo::apply(const Point &)
    //! Add ox and oy.
    GtransfoLinShift(double ox = 0., double oy = 0.) : GtransfoLin(ox, oy, 1., 0., 0., 1.) {}
    GtransfoLinShift(const Point &point) : GtransfoLin(point.x, point.y, 1., 0., 0., 1.){};
    double fit(const StarMatchList &starMatchList);

    int getNpar() const { return 2; }
};

/*=============================================================*/
//! just here to provide a specialized constructor, and fit.
class GtransfoLinRot : public GtransfoLin {
public:
    using Gtransfo::apply;  // to unhide apply(const Point&)

    GtransfoLinRot() : GtransfoLin(){};
    GtransfoLinRot(const double angleRad, const Point *center = nullptr, const double scaleFactor = 1.0);
    double fit(const StarMatchList &starMatchList);

    int getNpar() const { return 4; }
};

/*=============================================================*/

//! just here to provide specialized constructors. GtransfoLin fit routine.
class GtransfoLinScale : public GtransfoLin {
public:
    using Gtransfo::apply;  // to unhide apply(const Point&)
    //!
    GtransfoLinScale(const double scale = 1) : GtransfoLin(0.0, 0.0, scale, 0., 0., scale){};
    //!
    GtransfoLinScale(const double scaleX, const double scaleY)
            : GtransfoLin(0.0, 0.0, scaleX, 0., 0., scaleY){};

    int getNpar() const { return 2; }
};

/*==================WCS's transfo's =====================================*/

class BaseTanWcs : public Gtransfo {
protected:
    GtransfoLin linPix2Tan;  // pixels to tangent plane (internally in radians)
    std::unique_ptr<GtransfoPoly> corr;
    double ra0, dec0;   // in radians
    double cos0, sin0;  // cos(dec0), sin(dec0)

public:
    using Gtransfo::apply;  // to unhide apply(const Point&)

    BaseTanWcs(const GtransfoLin &pix2Tan, const Point &tangentPoint,
               const GtransfoPoly *corrections = nullptr);

    BaseTanWcs(const BaseTanWcs &original);

    void operator=(const BaseTanWcs &original);

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! The tangent point (in degrees)
    Point getTangentPoint() const;

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin getLinPart() const;

    //! the "correction" (non-owning pointer)
    const GtransfoPoly *getCorr() const { return corr.get(); }

    //! Assign the correction polynomial (what it means is left to derived classes)
    void setCorrections(std::unique_ptr<GtransfoPoly> corrections);

    //! the CRPIX values (this is WCS jargon), in 0-based coordinates
    Point getCrPix() const;

    //! transfo from pix to tangent plane (defined by derived classes)
    virtual GtransfoPoly getPix2TangentPlane() const = 0;

    //! Transforms from pixel space to tangent plane. deferred to actual implementations
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const = 0;

    ~BaseTanWcs();
};

class TanRaDec2Pix;  // the inverse of TanPix2RaDec.

//! the transformation that handles pix to sideral transfos (Gnomonic, possibly with polynomial distortions).
class TanPix2RaDec : public BaseTanWcs {
public:
    using Gtransfo::apply;  // to unhide apply(const Point&)
    //! pix2Tan describes the transfo from pix to tangent plane (in degrees). TangentPoint in degrees.
    //! Corrections are applied between Lin and deprojection parts (as in Swarp).
    TanPix2RaDec(const GtransfoLin &pix2Tan, const Point &tangentPoint,
                 const GtransfoPoly *corrections = nullptr);

    //! the transformation from pixels to tangent plane (coordinates in degrees)
    GtransfoPoly getPix2TangentPlane() const;

    //! transforms from pixel space to tangent plane
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const;

    TanPix2RaDec();

    //! composition with GtransfoLin
    TanPix2RaDec operator*(const GtransfoLin &right) const;

    std::unique_ptr<Gtransfo> reduceCompo(const Gtransfo *right) const;

    //! approximate inverse : it ignores corrections;
    TanRaDec2Pix invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    std::unique_ptr<Gtransfo> roughInverse(const Frame &region) const;

    //! Inverse transfo: returns a TanRaDec2Pix if there are no corrections, or the iterative solver if there
    //! are.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, const Frame &region) const;

    std::unique_ptr<Gtransfo> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(const StarMatchList &starMatchList);
};

//! Implements the (forward) SIP distorsion scheme
class TanSipPix2RaDec : public BaseTanWcs {
public:
    //! pix2Tan describes the transfo from pix to tangent plane (in degrees). TangentPoint in degrees.
    //! Corrections are applied before Lin.
    TanSipPix2RaDec(const GtransfoLin &pix2Tan, const Point &tangentPoint,
                    const GtransfoPoly *corrections = nullptr);

    //! the transformation from pixels to tangent plane (coordinates in degrees)
    GtransfoPoly getPix2TangentPlane() const;

    //! transforms from pixel space to tangent plane
    virtual void pix2TP(double xPixel, double yPixel, double &xTangentPlane, double &yTangentPlane) const;

    TanSipPix2RaDec();

    //! Inverse transfo: returns a TanRaDec2Pix if there are no corrections, or the iterative solver if there
    //! are.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, const Frame &region) const;

    std::unique_ptr<Gtransfo> clone() const;

    void dump(std::ostream &stream) const;

    //! Not implemented yet, because we do it otherwise.
    double fit(const StarMatchList &starMatchList);
};

//! This one is the Tangent Plane (called gnomonic) projection (from celestial sphere to tangent plane)
/*! this transfo does not implement corrections, since
   they are defined the other way around (from pixels to sky),
   and not invertible analytically. The inversion of tangent
   point WCS (TanPix2RaDec) is obtained via inverseTransfo().
*/

class TanRaDec2Pix : public Gtransfo {
    double ra0, dec0;  // tangent point (internally in radians)
    double cos0, sin0;
    GtransfoLin linTan2Pix;  // tangent plane to pixels (internally in radians)

public:
    using Gtransfo::apply;  // to unhide apply(const Point&)

    //! assume degrees everywhere.
    TanRaDec2Pix(const GtransfoLin &tan2Pix, const Point &tangentPoint);

    //!
    TanRaDec2Pix();

    //! The Linear part (corresponding to CD's and CRPIX's)
    GtransfoLin getLinPart() const;

    //! Resets the projection (or tangent) point
    void setTangentPoint(const Point &tangentPoint);

    //! tangent point coordinates (in degrees)
    Point getTangentPoint() const;

    //!
    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    //! transform with analytical derivatives
    void transformPosAndErrors(const FatPoint &in, FatPoint &out) const;

    //! exact typed inverse:
    TanPix2RaDec invert() const;

    //! Overload the "generic routine" (available for all Gtransfo types
    std::unique_ptr<Gtransfo> roughInverse(const Frame &region) const;

    //! Inverse transfo: returns a TanPix2RaDec.
    std::unique_ptr<Gtransfo> inverseTransfo(const double precision, const Frame &region) const;

    void dump(std::ostream &stream) const;

    std::unique_ptr<Gtransfo> clone() const;

    double fit(const StarMatchList &starMatchList);
};

//! signature of the user-provided routine that actually does the coordinate transfo for UserTransfo.
typedef void(GtransfoFun)(const double, const double, double &, double &, const void *);

//! a run-time transfo that allows users to define a Gtransfo with minimal coding (just the transfo routine).
class UserTransfo : public Gtransfo {
private:
    GtransfoFun *_userFun;
    const void *_userData;

public:
    using Gtransfo::apply;  // to unhide apply(const Point&)

    //! the transfo routine and extra data that it may need.
    UserTransfo(GtransfoFun &fun, const void *userData);

    void apply(const double xIn, const double yIn, double &xOut, double &yOut) const;

    void dump(std::ostream &stream = std::cout) const;

    double fit(const StarMatchList &starMatchList);

    std::unique_ptr<Gtransfo> clone() const;
};

//! The virtual constructor from a file
std::unique_ptr<Gtransfo> gtransfoRead(const std::string &fileName);
//! The virtual constructor from a file
std::unique_ptr<Gtransfo> gtransfoRead(std::istream &s);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_GTRANSFO_H
