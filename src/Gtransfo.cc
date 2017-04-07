#include <iostream>
#include <iomanip>
#include <iterator> /* for ostream_iterator */
#include <math.h> // for sin and cos and may be others
#include <fstream>
#include "assert.h"
#include <sstream>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/pex/exceptions.h"
#include "Eigen/Cholesky"

namespace pexExcept = lsst::pex::exceptions;

using namespace std;

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.Gtransfo");
}

namespace lsst {
namespace jointcal {

bool IsIdentity(const Gtransfo *a_transfo)
{ return (dynamic_cast<const GtransfoIdentity*>(a_transfo) != nullptr);}

static double sqr(double x) {return x*x;}

bool IsIntegerShift(const Gtransfo *a_transfo)
{
  const GtransfoPoly* shift = dynamic_cast<const GtransfoPoly*>(a_transfo);
  if (shift == nullptr) return false;

  static const double eps = 1e-5;

  double dx = shift->Coeff(0,0,0);
  double dy = shift->Coeff(0,0,1);

  static Point dumb(4000,4000);
  if (fabs(dx - int(floor(dx+0.5))) < eps &&
      fabs(dy - int(floor(dy+0.5))) < eps &&
      fabs(dumb.x+dx - shift->apply(dumb).x) < eps &&
      fabs(dumb.y+dy - shift->apply(dumb).y) < eps)
    return true;

  return false;
}

/********* Gtransfo ***********************/

std::unique_ptr<Gtransfo> Gtransfo::ReduceCompo(const Gtransfo *Right) const
{// by default no way to compose
  if (Right) {} // avoid a warning
  return std::unique_ptr<Gtransfo>(nullptr);
}


double Gtransfo::Jacobian(const double X, const double Y) const
{
double x2,y2;
double eps=X*0.01;
if (eps == 0) eps = 0.01;
apply(X,Y, x2, y2);
double dxdx, dydx;
apply(X+eps, Y, dxdx, dydx);
dxdx -= x2; dydx -= y2;
double dxdy, dydy;
apply(X, Y+eps, dxdy, dydy);
dxdy -= x2; dydy -= y2;
return ((dxdx * dydy - dxdy * dydx)/(eps*eps));
}

/*! the Derivative is represented by a GtransfoLin, in which
  (hopefully), the offset terms are zero. Derivative should
  transform a vector of offsets into a vector of offsets. */
void Gtransfo::Derivative(const Point &Where,
			  GtransfoLin &Der, const double Step) const
{
  double x = Where.x;
  double y = Where.y;
  double xp0,yp0;
  apply(x, y, xp0, yp0);

  double xp,yp;
  apply(x+Step, y, xp, yp);
  Der.a11() = (xp-xp0)/Step;
  Der.a21() = (yp-yp0)/Step;
  apply(x, y + Step, xp, yp);
  Der.a12() = (xp-xp0)/Step;
  Der.a22() = (yp-yp0)/Step;
  Der.dx() = 0;
  Der.dy() = 0;
}

GtransfoLin Gtransfo::LinearApproximation(const Point &Where,
					  const double Step) const
{
  Point outWhere = apply(Where);
  GtransfoLin der;
  Derivative(Where, der, Step);
  return GtransfoLinShift(outWhere.x, outWhere.y)*der*GtransfoLinShift(-Where.x, -Where.y);
}


void Gtransfo::TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const
{
  FatPoint res; // in case In and Out are the same address...
  res = apply(In);
  GtransfoLin der;
  // could save a call here, since Derivative needs the transform of where that we already have
  // 0.01 may not be a very good idea in all cases. May be we should provide a way of altering that.
  Derivative(In, der, 0.01);
  double a11 = der.A11();
  double a22 = der.A22();
  double a21 = der.A21();
  double a12 = der.A12();
  res.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  res.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  res.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;
  Out = res;
}


void Gtransfo::TransformErrors(const Point &Where,
			       const double *VIn, double *VOut) const
{
  GtransfoLin der;
  Derivative(Where, der, 0.01);
  double a11 = der.A11();
  double a22 = der.A22();
  double a21 = der.A21();
  double a12 = der.A12();

  /*      (a11 a12)          (vxx  vxy)
     M =  (       )  and V = (        )
          (a21 a22)          (xvy  vyy)

     Vxx = Vin[0], vyy = Vin[1], Vxy = Vin[2];
     we want to compute M*V*tp(M)
     A lin alg light package would be perfect...
  */
  int xx = 0;
  int yy = 1;
  int xy = 2;
  // M*V :

  double b11 = a11 * VIn[xx] + a12 * VIn[xy];
  double b22 = a21 * VIn[xy] + a22 * VIn[yy];
  double b12 = a11 * VIn[xy] + a12 * VIn[yy];
  double b21 = a21 * VIn[xx] + a22 * VIn[xy];

  // (M*V) * tp(M)

  VOut[xx] = b11 * a11 + b12 * a12;
  VOut[xy] = b11 * a21 + b12 * a22;
  VOut[yy] = b21 * a21 + b22 * a22;
}

std::unique_ptr<Gtransfo> Gtransfo::RoughInverse(const Frame &Region) const
{
  // "in" and "out" refer to the inverse direction.
  Point centerOut = Region.Center();
  Point centerIn = apply(centerOut);
  GtransfoLin der;
  Derivative(centerOut,der,sqrt(Region.Area())/5.);
  der = der.invert();
  der = GtransfoLinShift(centerOut.x, centerOut.y)
    *der
    *GtransfoLinShift(-centerIn.x, -centerIn.y);
  return std::unique_ptr<Gtransfo>(new GtransfoLin(der));
}


/* implement one in Gtransfo, so that all derived
   classes do not need to provide one... */


/* the routines that follow are used for ea generic parameter
   transformation serialization, used e.g. for fits. Enables
   to manipulate transformation parameters as vectors.
*/


// not dummy : what it does is virtual because ParamRef is virtual.
void Gtransfo::GetParams(double *Params) const
{
  int npar = Npar();
  for (int i=0; i<npar ; ++i) Params[i] = ParamRef(i);
}


void Gtransfo::OffsetParams(const double *Params)
{
  int npar = Npar();
  for (int i=0; i<npar ; ++i) ParamRef(i) += Params[i];
}

double Gtransfo::ParamRef(const int) const
{
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, std::string("Gtransfo::ParamRef should never be called "));
}

double &Gtransfo::ParamRef(const int)
{
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Gtransfo::ParamRef should never be called ");
}

void Gtransfo::ParamDerivatives(const Point &Where,
				double *Dx, double *Dy) const
{
  if ((Where.x || Dx ) || Dy) {} // compilation warning killer
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Gtransfo::ParamDerivatives() should never be called ");
}


ostream & operator << (ostream &stream, const Gtransfo & T)
           {T.dump(stream); return stream;}



void Gtransfo::Write(const std::string &FileName) const
{
  ofstream s(FileName.c_str());
  Write(s);
  bool ok = !s.fail();
  s.close();
  if (!ok)
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Gtransfo::Write, something went wrong for file " + FileName );
}

void Gtransfo::Write(ostream &stream) const
{
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "Gtransfo::Write(ostream), should never be called. MEans that it is missing in some derived class ");
}


/******************* GTransfoInverse ****************/
/* inverse transformation, solved by iterations. Before using
   it (probably via Gtransfo::InverseTransfo), consider
   seriously StarMatchList::InverseTransfo */
class GtransfoInverse : public Gtransfo {

private:
  std::unique_ptr<Gtransfo> direct;
  std::unique_ptr<Gtransfo> roughInverse;
  double precision2;


public:
  GtransfoInverse(const Gtransfo* Direct,
		  const double Precision,
		  const Frame& Region);

  //! implements an iterative (Gauss-Newton) solver. It resorts to the Derivative function: 4 calls to the direct transfo per iteration.
  void apply(const double Xin, const double Yin,
	     double &Xout, double  &Yout) const;

  void dump(ostream &stream) const;

  double fit(const StarMatchList &List);

  virtual std::unique_ptr<Gtransfo> Clone() const;

  GtransfoInverse(const GtransfoInverse&);

  //! Overload the "generic routine"
  std::unique_ptr<Gtransfo> RoughInverse(const Frame &) const
  {
    return direct->Clone();
  }

  //! Inverse transfo: returns the direct one!
  std::unique_ptr<Gtransfo> InverseTransfo(const double,
			   const Frame &) const
  {
    return direct->Clone();
  }

  ~GtransfoInverse();

private:
  void operator = (const GtransfoInverse &);

};

std::unique_ptr<Gtransfo> Gtransfo::InverseTransfo(const double Precision,
				   const Frame& Region) const
{
  return std::unique_ptr<Gtransfo>(new GtransfoInverse(this,Precision,Region));
}


GtransfoInverse::GtransfoInverse(const Gtransfo* Direct,
				 const double Precision,
				 const Frame& Region)

{
  direct = Direct->Clone();
  roughInverse = Direct->RoughInverse(Region);
  precision2 = Precision*Precision;
}

GtransfoInverse::GtransfoInverse(const GtransfoInverse& Model) : Gtransfo()
{
  direct = Model.direct->Clone();
  roughInverse = Model.roughInverse->Clone();
  precision2 = Model.precision2;
}

GtransfoInverse::~GtransfoInverse() { }

void GtransfoInverse::operator = (const GtransfoInverse & Model)
{
  direct = Model.direct->Clone();
  roughInverse = Model.roughInverse->Clone();
  precision2 = Model.precision2;
}

void GtransfoInverse::apply(const double Xin, const double Yin,
			    double &Xout, double  &Yout) const
{
  Point in(Xin,Yin);
  Point outGuess = roughInverse->apply(in);
  GtransfoLin directDer, reverseDer;
  int loop = 0;
  int maxloop = 20;
  double move2;
  do
    {
      loop++;
      Point inGuess = direct->apply(outGuess);
      direct->Derivative(outGuess,directDer);
      reverseDer = directDer.invert();
      double xShift, yShift;
      reverseDer.apply(Xin - inGuess.x, Yin - inGuess.y, xShift, yShift);
      outGuess.x += xShift;
      outGuess.y += yShift;
      move2 = xShift*xShift+yShift*yShift;
    } while (( move2 > precision2) && (loop < maxloop));
  if (loop == maxloop)
    LOGLS_WARN(_log, "Problems applying GtransfoInverse at " << in);
  Xout = outGuess.x;
  Yout = outGuess.y;
}

void GtransfoInverse::dump(ostream &stream) const
{
  stream << " GtransfoInverse of  :" << endl
	 << *direct << endl;
}

double GtransfoInverse::fit(const StarMatchList &)
{
  throw pexExcept::RuntimeError("Cannot fit a GtransfoInverse. Use StarMatchList::inverseTransfo instead.");
}

std::unique_ptr<Gtransfo> GtransfoInverse::Clone() const
{
  return std::unique_ptr<Gtransfo>(new GtransfoInverse(*this));
}


/************* GtransfoComposition **************/


// This class was done to allow composition of Gtransfo's, without specifications of their types.
// does not need to be public. Invoked  by GtransfoCompose(Left,Right)



//! Private class to handle Gtransfo compositions (i.e. piping). Use the routine GtransfoCompose if you need this functionnality.
class GtransfoComposition : public Gtransfo {
  private :
      std::unique_ptr<Gtransfo> first, second;
  public :
    //! will pipe transfos
    GtransfoComposition(const Gtransfo *Second, const Gtransfo *First);

    //! return Second(First(Xin,Yin))
    void apply(const double Xin, const double Yin, double &Xout, double &Yout) const;
    void dump(ostream &stream = cout) const;

    //!
    double fit(const StarMatchList &List);

    std::unique_ptr<Gtransfo> Clone() const;
    ~GtransfoComposition();
};


GtransfoComposition::GtransfoComposition(const Gtransfo *Second, const Gtransfo *First)
{
first =  First->Clone();
second = Second->Clone();
}

void GtransfoComposition::apply(const double Xin, const double Yin, double &Xout, double &Yout) const
{
double xout,yout;
first->apply(Xin,Yin, xout,yout);
second->apply(xout,yout,Xout,Yout);
}

void GtransfoComposition::dump(ostream &stream) const
{
first->dump(stream); second->dump(stream);
}

double GtransfoComposition::fit(const StarMatchList &List)
{
  /* fits only one of them. could check that first can actually be fitted... */
  return first->fit(List);
}

std::unique_ptr<Gtransfo> GtransfoComposition::Clone() const
{
return std::unique_ptr<Gtransfo>(new GtransfoComposition(second.get(),first.get()));
}

GtransfoComposition::~GtransfoComposition() { }

/*!  This routine implements "run-time" compositions. When
 there is a possible "reduction" (e.g. compositions of polynomials),
 GtransfoCompose detects it and returns a genuine Gtransfo.
 */
std::unique_ptr<Gtransfo> GtransfoCompose(const Gtransfo *Left, const Gtransfo *Right)
{
  /* is Right Identity ? if Left is Identity , GtransfoIdentity::ReduceCompo does the right job */
  if (IsIdentity(Right))
    {
      return Left->Clone();
    }
  /* Try to use the ReduceCompo method from Left. If absent,
     Gtransfo::ReduceCompo return NULL. ReduceCompo is non trivial for
     polynomials */
  std::unique_ptr<Gtransfo> composition(Left->ReduceCompo(Right));
  /* composition == NULL means no reduction : just build a Composition
     that pipelines "Left" and "Right" */
  if (composition == nullptr) return std::unique_ptr<Gtransfo>(new GtransfoComposition(Left,Right));
  else return composition;
}


// just a speed up, to avoid useless numerical derivation.
void GtransfoIdentity::Derivative(const Point &,
				  GtransfoLin &Derivative,
				  const double) const
{
  Derivative = GtransfoLin();
}


GtransfoLin GtransfoIdentity::LinearApproximation(const Point &, const double) const
{
  GtransfoLin result;
  return result; // rely on default Gtransfolin constructor;
}


void GtransfoIdentity::Write(ostream &s) const
{
  s << "GtransfoIdentity 1" << endl;
}




void GtransfoIdentity::Read(istream &s)
{
  int format;
  s >> format;
  if (format != 1)
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, " GtransfoIdentity::Read : format is not 1 " );
}


/***************  GtransfoPoly **************************************/



//! Default transfo : identity for all degrees (>=1 )

GtransfoPoly::GtransfoPoly(const unsigned Deg)
{
  deg = Deg;
  nterms = (deg+1)*(deg+2)/2;

  // allocate and fill coefficients
  coeffs.resize(2*nterms,0.);
  // the default is supposed to be the identity, (for deg>=1).
  if (deg>=1)
    {
      Coeff(1,0,0) = 1;
      Coeff(0,1,1) = 1;
    }
}


  //#ifdef TO_BE_FIXED
GtransfoPoly::GtransfoPoly(const Gtransfo* T,
			   const Frame& F,
			   unsigned Degree,
			   unsigned NPoint)
{
  StarMatchList sm;

  double step = sqrt(fabs(F.Area())/double(NPoint));
  for (double x=F.xMin+step/2; x<=F.xMax; x+=step)
  for (double y=F.yMin+step/2; y<=F.yMax; y+=step)
    {
      auto pix = std::make_shared<BaseStar>(x,y,0);
      double xtr, ytr;
      T->apply(x,y,xtr,ytr);
      auto tp = std::make_shared<BaseStar>(xtr,ytr,0);
      /* These are fake stars so no need to transform fake errors.
	 all errors (and weights) will be equal : */
      sm.push_back(StarMatch(*pix,*tp,pix,tp));
    }
  GtransfoPoly ret(Degree);
  ret.fit(sm);
  *this = ret;
}
  //#endif






void GtransfoPoly::compute_monomials(double Xin, double Yin,
				     double *Monom) const
{
  /* The ordering of monomials is implemented here.
     You may not change it without updating the "mapping" routines
    Coeff(unsigned, unsigned, unsigned).
    I (P.A.) did not find a clever way to loop over monomials.
    Improvements welcome.
    This routine is used also by the fit to fill monomials.
    We could certainly be more elegant.
  */

  double xx = 1;
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      double yy = 1;
      unsigned k=ix*(ix+1)/2;
      for (unsigned iy = 0; iy<=deg-ix; ++iy)
	{
	  //	  assert(k<nterms);
	Monom[k] = xx*yy;
	yy *= Yin;
	k+= ix+iy+2;
      }
    xx *= Xin;
    }
}


void GtransfoPoly::SetDegree(const unsigned Deg)
{
  deg = Deg;
  unsigned old_nterms = nterms;
  nterms = (deg+1)*(deg+2)/2;

  // temporarily save coefficients
  vector<double> old_coeffs = coeffs;
  // reallocate enough size
  coeffs.resize(2*nterms);
  // reassign to zero (this is necessary because ycoeffs
  // are after xcoeffs and so their meaning changes
  for (unsigned k=0; k<nterms; ++k) coeffs[k] = 0;
  // put back what we had before
  unsigned kmax = min(old_nterms, nterms);
  for (unsigned k=0; k<kmax; ++k)
    {
      coeffs[k] = old_coeffs[k];  // x terms
      coeffs[k+nterms] = old_coeffs[k+old_nterms]; // y terms
    }
}

/* this is reasonably fast, when optimized */
void GtransfoPoly::apply(const double Xin, const double Yin,
			 double &Xout, double &Yout) const
{

  /*
    This routine computes the monomials only once for both
    polynomials.  This is why GtransfoPoly does not use an auxilary
    class (such as PolyXY) to handle each polynomial.

    The code works even if &Xin == &Xout (or &Yin == &Yout)
    It uses Variable Length Allocation (VLA) rather than a vector<double>
    because allocating the later costs about 50 ns. All VLA uses are tagged.
  */
  double monomials[nterms];    // this is VLA, which is (perhaps) not casher C++
  compute_monomials(Xin, Yin, monomials);

  Xout = 0;
  Yout = 0;
  const double *c = &coeffs[0];
  const double *pm = &monomials[0];
  // the ordering of the coefficients and the monomials are identical.
  for (int k=nterms; k--; ) Xout +=  (*(pm++))*(*(c++));
  pm = &monomials[0];
  for (int k=nterms; k--; ) Yout +=  (*(pm++))*(*(c++));
}


void GtransfoPoly::Derivative(const Point &Where,
			      GtransfoLin &Der, const double Step) const
{ /* routine checked against numerical derivatives from Gtransfo::Derivative */
  if (deg == 1)
    {
      Der = GtransfoLin(*this);
      Der.dx() = Der.dy() = 0;
      return;
    }

  double dermx[2*nterms];   //VLA
  double *dermy = dermx+nterms;
  double xin = Where.x;
  double yin = Where.y;

  double xx = 1;
  double xxm1 = 1; // xx^(ix-1)
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      unsigned k=(ix)*(ix+1)/2;
      // iy = 0
      dermx[k] = ix*xxm1;
      dermy[k] = 0;
      k+= ix+2;
      double yym1 = 1; // yy^(iy-1)
      for (unsigned iy = 1; iy<=deg-ix; ++iy)
	{
	  dermx[k] = ix*xxm1*yym1*yin;
	  dermy[k] = iy*xx*yym1;
	  yym1 *= yin;
	  k+= ix+iy+2;
	}
    xx *= xin;
    if (ix>=1) xxm1 *= xin;
    }

  Der.dx() = 0;
  Der.dy() = 0;

  const double *mx = &dermx[0];
  const double *my = &dermy[0];
  const double *c = &coeffs[0];
  // dx'
  double a11=0, a12 = 0;
  for (int k=nterms; k--; )
    {
      a11 += (*(mx++))*(*c);
      a12 += (*(my++))*(*(c++));
    }
  Der.a11() = a11;
  Der.a12() = a12;
  // dy'
  double a21 = 0, a22 = 0;
  mx = &dermx[0];
  my = &dermy[0];
  for (int k=nterms; k--; )
    {
      a21 += (*(mx++))*(*c);
      a22 += (*(my++))*(*(c++));
    }
  Der.a21() = a21;
  Der.a22() = a22;

}

void GtransfoPoly::TransformPosAndErrors(const FatPoint &In, FatPoint &Out) const
{
  /*
     The results from this routine were compared to what comes out
     from apply and TransformErrors. The Derivative routine was
     checked against numerical derivatives from
     Gtransfo::Derivative. (P.A dec 2009).

     This routine could be made much simpler by calling apply and
     Derivative (i.e. you just suppress it, and the fallback is the
     generic version in Gtransfo).  BTW, I checked that both routines
     provide the same result. This version is however faster
     (monomials get recycled).
  */
  double monomials[nterms]; // VLA

  FatPoint  res; // to store the result, because nothing forbids &In == &Out.

  double dermx[2*nterms]; // monomials for derivative w.r.t. x (VLA)
  double *dermy = dermx+nterms;  // same for y
  double xin = In.x;
  double yin = In.y;

  double xx = 1;
  double xxm1 = 1; // xx^(ix-1)
  for (unsigned ix = 0; ix<=deg; ++ix)
    {
      unsigned k=(ix)*(ix+1)/2;
      // iy = 0
      dermx[k] = ix*xxm1;
      dermy[k] = 0;
      monomials[k] = xx;
      k+= ix+2;
      double yy = yin;
      double yym1 = 1; // yy^(iy-1)
      for (unsigned iy = 1; iy<=deg-ix; ++iy)
	{
	  monomials[k] = xx*yy;
	  dermx[k] = ix*xxm1*yy;
	  dermy[k] = iy*xx*yym1;
	  yym1 *= yin;
	  yy *= yin;
	  k+= ix+iy+2;
	}
    xx *= xin;
    if (ix>=1) xxm1 *= xin;
    }

  // output position
  double xout = 0, yout=0;
  const double *c = &coeffs[0];
  const double *pm = &monomials[0];
  for (int k=nterms; k--; ) xout +=  (*(pm++))*(*(c++));
  pm = &monomials[0];
  for (int k=nterms; k--; ) yout +=  (*(pm++))*(*(c++));
  res.x = xout; res.y = yout;

  // derivatives
  c = &coeffs[0];
  const double *mx = &dermx[0];
  const double *my = &dermy[0];
  double a11=0, a12 = 0;
  for (int k=nterms; k--; )
    {
      a11 += (*(mx++))*(*c);
      a12 += (*(my++))*(*(c++));
    }

  double a21 = 0, a22 = 0;
  mx = &dermx[0];
  my = &dermy[0];
  for (int k=nterms; k--; )
    {
      a21 += (*(mx++))*(*c);
      a22 += (*(my++))*(*(c++));
    }

  // output co-variance
  res.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  res.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  res.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;
  Out = res;
}


/* The coefficient ordering is defined both here *AND* in the
   GtransfoPoly::apply, GtransfoPoly::Derivative, ... routines
   Change all or none ! */

double GtransfoPoly::Coeff(const unsigned Degx, const unsigned Degy,
			   const unsigned WhichCoord) const
{
  assert((Degx+Degy<=deg) && WhichCoord<2);
  /* this assertion above is enough to ensure that the index used just
     below is within bounds since the reserved length is
     2*nterms=(deg+1)*(deg+2) */
  return coeffs[(Degx+Degy)*(Degx+Degy+1)/2+Degy+WhichCoord*nterms];
}


double& GtransfoPoly::Coeff(const unsigned Degx, const unsigned Degy,
			   const unsigned WhichCoord)
{
  assert((Degx+Degy<=deg) && WhichCoord<2);
  return coeffs[(Degx+Degy)*(Degx+Degy+1)/2+Degy+WhichCoord*nterms];
}

double GtransfoPoly::CoeffOrZero(const unsigned Degx, const unsigned Degy,
				 const unsigned WhichCoord) const
{
  //  assert((Degx+Degy<=deg) && WhichCoord<2);
  assert(WhichCoord<2);
  if (Degx+Degy<=deg)
    return coeffs[(Degx+Degy)*(Degx+Degy+1)/2+Degy+WhichCoord*nterms];
  return 0;
}

/* parameter serialization for "virtual" fits */
double GtransfoPoly::ParamRef(const int i) const
{
  assert(unsigned(i)<2*nterms);
  return coeffs[i];
}


double& GtransfoPoly::ParamRef(const int i)
{
  assert(unsigned(i)<2*nterms);
  return coeffs[i];
}

void GtransfoPoly::ParamDerivatives(const Point &Where,
				    double *Dx, double *Dy) const
{/* first half : dxout/dpar, second half : dyout/dpar */
  compute_monomials(Where.x, Where.y, Dx );
  for (unsigned k=0; k<nterms; ++k)
    {
      Dy[nterms+k] = Dx[k];
      Dx[nterms+k] = Dy[k] = 0;
    }
}

/*
  mapping coefficients with names (inherited from GtransfoLin,GtransfoQuad,GtransfoCub} */
typedef struct
{
  const char *name;
  const unsigned char px,py,whichCoord;
} CoeffTagStruct;

static const CoeffTagStruct CoeffTags[] =
  {
    {"dx",   0,0,0},
    {"dy",   0,0,1},
    {"a11",  1,0,0},
    {"a12",  0,1,0},
    {"a22",  0,1,1},
    {"a21",  1,0,1},
    {"a1x2", 2,0,0},
    {"a1xy", 1,1,0},
    {"a1y2", 0,2,0},
    {"a2y2", 0,2,1},
    {"a2xy", 1,1,1},
    {"a2x2", 2,0,1},
    {"a1x3", 3,0,0},
    {"a1x2y",2,1,0},
    {"a1xy2",1,2,0},
    {"a1y3", 0,3,0},
    {"a2y3", 0,3,1},
    {"a2xy2",1,2,1},
    {"a2x2y",2,1,1},
    {"a2x3", 3,0,1}
  };

static unsigned NCoeffTags = sizeof(CoeffTags)/sizeof(CoeffTags[0]);



static unsigned tag_pos(const char *Name)
{
  for (unsigned k=0; k<NCoeffTags; ++ k)
    {
      const CoeffTagStruct &t=CoeffTags[k];
      if (strcmp(t.name, Name) == 0)
	return k;
    }
  //  stringstream message;
  //  message << "GtransfoPoly::Coeff(const char *Name) : unknown name : \""
  //	  << string(Name) << '\"';
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
		    "GtransfoPoly::Coeff(const char *Name) : unknown name : \""+std::string(Name)+"\"");
}

bool GtransfoPoly::HasCoeff(const char* Name) const
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];
  return (t.px+t.py<=deg);
}

double GtransfoPoly::Coeff(const char *Name) const
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];
  return Coeff(t.px, t.py, t.whichCoord);
}

double& GtransfoPoly::Coeff(const char *Name)
{
  const CoeffTagStruct &t=CoeffTags[tag_pos(Name)];
  return Coeff(t.px, t.py, t.whichCoord);
}

/* utility for the dump(ostream&) routine */

static string monomial_string(const unsigned powx, const unsigned powy)
{
  stringstream ss;
  if (powx+powy) ss<<"*";
  if (powx>0) ss<< "x";
  if (powx>1) ss<<"^"<<powx;
  if (powy>0) ss<< "y";
  if (powy>1) ss<<"^"<<powy;
  return ss.str();
}

void  GtransfoPoly::dump(ostream &S) const
{
  for (unsigned ic=0; ic<2; ++ic)
    {
      if (ic==0)   S << "newx = ";
      else S << "newy = ";
      for (unsigned p = 0; p<=deg; ++p)
	for (unsigned py=0; py<=p; ++py)
	  {
	    if (p+py != 0) S<< " + ";
	    S << Coeff(p-py,py,ic) << monomial_string(p-py,py);
	  }
      S << endl;
    }
  if (deg>0)
    S << " Linear Determinant = " << Determinant() << endl ;
}

double GtransfoPoly::Determinant() const
{
  if (deg >= 1)
    return Coeff(1,0,0)*Coeff(0,1,1) - Coeff(0,1,0)*Coeff(1,0,1);
  return 0;
}

//! Returns the transformation that maps the input frame along both axes to [-1,1]
GtransfoLin NormalizeCoordinatesTransfo(const Frame & F)
{
  Point center = F.Center();
  return GtransfoLinScale(2./F.Width(), 2./F.Height())*GtransfoLinShift(-center.x,-center.y);
}



/*utility for the GtransfoPoly::fit() routine */
static GtransfoLin shift_and_normalize(const StarMatchList &List)
{
  double xav=0;
  double x2=0;
  double yav=0;
  double y2=0;
  double count=0;
  for (auto it = List.begin(); it != List.end(); ++it)
    {
      const StarMatch &a_match = *it;
      const Point &point1 = a_match.point1;
      xav += point1.x;
      yav += point1.y;
      x2 += sqr(point1.x);
      y2 += sqr(point1.y);
      count++;
    }
  if (count==0) return GtransfoLin();
  xav /= count;
  yav /= count;
  // 3.5 stands for sqrt(12).
  double xspan = 3.5*sqrt(x2/count-sqr(xav));
  double yspan = 3.5*sqrt(y2/count-sqr(yav));
  return GtransfoLinScale(2./xspan, 2./yspan) * GtransfoLinShift(-xav,-yav);
}

static double sq(double x) { return x*x;}

double GtransfoPoly::do_the_fit(const StarMatchList &List,
				const Gtransfo &ShiftToCenter,
				const bool UseErrors)
{
  Eigen::MatrixXd A(2*nterms,2*nterms);   A.setZero();
  Eigen::VectorXd B(2*nterms);   B.setZero();
  double sumr2 = 0;
  double monomials[nterms];
  for (auto it = List.begin(); it != List.end(); ++it)
    {
      const StarMatch &a_match = *it;
      Point tmp = ShiftToCenter.apply(a_match.point1);
      FatPoint point1(tmp, a_match.point1.vx, a_match.point1.vy,
		      a_match.point1.vxy);
      const FatPoint &point2 = a_match.point2;
      double wxx,wyy,wxy;
      FatPoint tr1;
      compute_monomials(point1.x, point1.y, monomials);
      if (UseErrors)
	{
	  TransformPosAndErrors(point1, tr1);// we might consider recycling the monomials
	  double vxx = (tr1.vx+point2.vx);
	  double vyy = (tr1.vy+point2.vy);
	  double vxy = (tr1.vxy+point2.vxy);
	  double det = vxx*vyy-vxy*vxy;
	  wxx = vyy/det;
	  wyy = vxx/det;
	  wxy = -vxy/det;
	}
      else
	{
	  wxx = wyy = 1; wxy = 0;
	  apply(point1.x, point1.y ,tr1.x, tr1.y);
	}
      double resx = point2.x - tr1.x;
      double resy = point2.y - tr1.y;
      sumr2 += wxx*sq(resx) + wyy*sq(resy)
	+2*wxy*resx*resy;

      double bxcoeff = wxx*resx + wxy*resy ;
      double bycoeff = wyy*resy + wxy*resx;
      for (unsigned j=0; j<nterms; ++j)
	{
	  for (unsigned i=j; i<nterms; ++i)
	    {
	      A(i,j) += wxx*monomials[i]*monomials[j];
	      A(i+nterms,j+nterms) += wyy*monomials[i]*monomials[j];
	      A(j,i+nterms) = A(i,j+nterms) += wxy*monomials[i]*monomials[j];
	    }
	  B(j)        += bxcoeff*monomials[j];
	  B(j+nterms) += bycoeff*monomials[j];
	}
    } // end loop on points
  Eigen::LDLT<Eigen::MatrixXd, Eigen::Lower> factor(A);
  // should probably throw
  if (factor.info() != Eigen::Success)
    {
      LOGL_ERROR(_log, "GtransfoPoly::fit could not factorize");
      return -1;
    }

  Eigen::VectorXd sol = factor.solve(B);
  for (unsigned k=0; k< 2*nterms; ++k) coeffs[k] += sol(k);
  if (List.size() == nterms) return 0;
  return (sumr2-B.dot(sol));
}


double  GtransfoPoly::fit(const StarMatchList &List)
{
  if (List.size()< nterms)
    {
      LOGLS_FATAL(_log, "GtransfoPoly::fit trying to fit a polynomial transfo of degree " << deg << " with only " << List.size() << " matches.");
      return -1;
    }

  GtransfoPoly conditionner = shift_and_normalize(List);

  do_the_fit(List, conditionner, false); // get a rough solution
  do_the_fit(List, conditionner, true); // weight with it
  double chi2 = do_the_fit(List, conditionner, true); // once more

  (*this) = (*this)*conditionner;
  if (List.size() == nterms) return 0;
  return chi2;
}


std::unique_ptr<Gtransfo> GtransfoPoly::ReduceCompo(const Gtransfo *Right) const
{
  const GtransfoPoly *p = dynamic_cast<const GtransfoPoly *>(Right);
  if (p)
    {
      if (Degree() == 1 && p->Degree() == 1)
	return std::unique_ptr<Gtransfo>(new GtransfoLin((*this)*(*p))); // does the composition
      else
	return std::unique_ptr<Gtransfo>(new GtransfoPoly((*this)*(*p))); // does the composition
    }
  else return std::unique_ptr<Gtransfo>(nullptr);
}

/*  PolyXY the class used to perform polynomial algebra (and in
    particular composition) at the coefficient level. This class
    handles a single polynomial, while a GtransfoPoly is a couple of
    polynomials. This class does not have any routine to evaluate
    polynomials. Efficiency is not a concern since these routines are
    seldom used.  There is no need to expose this tool class to
    Gtransfo users.
*/


class PolyXY
{
  unsigned deg;
  unsigned nterms;
  vector<long double> coeffs;

public :

  PolyXY(const int Deg) : deg(Deg), nterms((deg+1)*(deg+2)/2)
  {
    coeffs.reserve(nterms);
    coeffs.insert(coeffs.begin(), nterms, 0L); // fill & initialize to 0.
  }

  unsigned Deg() const { return deg;}

  PolyXY(const GtransfoPoly &P, const unsigned WhichCoord)
    : deg(P.Degree()) , nterms((deg+1)*(deg+2)/2) , coeffs(nterms, 0L)
  {
    for (unsigned px=0; px<=deg; ++px)
      for (unsigned py=0; py<=deg-px; ++py)
	Coeff(px,py) = P.Coeff(px,py,WhichCoord);
  }


  long double Coeff(const unsigned powx, const unsigned powy) const
  {
    assert(powx+powy<=deg);
    return coeffs.at((powx+powy)*(powx+powy+1)/2+powy);
  }

  long double &Coeff(const unsigned powx, const unsigned powy)
  {
    assert(powx+powy<=deg);
    return coeffs.at((powx+powy)*(powx+powy+1)/2+powy);
  }

};


/* =====================  PolyXY Algebra routines ================== */

static void operator += (PolyXY &Left, const PolyXY &Right)
{
  unsigned rdeg = Right.Deg();
  assert(Left.Deg()>= rdeg);
  for (unsigned i=0; i<= rdeg; ++i)
    for (unsigned j=0; j<= rdeg-i; ++j)
      Left.Coeff(i,j) += Right.Coeff(i,j);
}

/* multiplication by a scalar */
static PolyXY operator * (const long double &a, const PolyXY &P)
{
  PolyXY result(P);
  // no direct access to coefficients: do it the soft way
  unsigned deg = P.Deg();
  for (unsigned i=0; i<=deg; ++i)
    for (unsigned j=0; j<=deg-i; ++j)
      result.Coeff(i,j) *= a;
  return result;
}


/*! result(x,y) = P1(x,y)*P2(x,y) */
static PolyXY Product(const PolyXY &P1, const PolyXY &P2)
{
  unsigned deg1 = P1.Deg();
  unsigned deg2 = P2.Deg();
  PolyXY result(deg1+deg2);
  for (unsigned i1=0; i1<=deg1; ++i1)
    for (unsigned j1=0; j1<=deg1-i1; ++j1)
      for (unsigned i2=0; i2<=deg2; ++i2)
	for (unsigned j2=0; j2<=deg2-i2; ++j2)
	  result.Coeff(i1+i2,j1+j2) += P1.Coeff(i1,j1)*P2.Coeff(i2,j2);
  return result;
}


/* Powers[k](x,y) = P(x,y)**k, 0 <= k <= MaxP */
static void ComputePowers(const PolyXY &P, const unsigned MaxP, vector<PolyXY> &Powers)
{
  Powers.reserve(MaxP+1);
  Powers.push_back(PolyXY(0)); Powers[0].Coeff(0,0) = 1L;
  for (unsigned k=1; k<=MaxP; ++k) Powers.push_back(Product(Powers[k-1],P));
}

/*! result(x,y) = P(Px(x,y),Py(x,y)) */
static PolyXY Composition(const PolyXY &P, const PolyXY &Px, const PolyXY &Py)
{
  unsigned pdeg = P.Deg();
  PolyXY result(pdeg*max(Px.Deg(), Py.Deg()));
  vector<PolyXY> PxPowers;
  vector<PolyXY> PyPowers;
  ComputePowers(Px, pdeg, PxPowers);
  ComputePowers(Py, pdeg, PyPowers);
  for (unsigned px=0 ; px <= pdeg; ++px)
    for (unsigned py=0; py <= pdeg-px; ++py)
      result += P.Coeff(px,py)*Product(PxPowers.at(px), PyPowers.at(py));
  return result;
}

/* ===================== end of  PolyXY Algebra routines ============= */


/* reducing polynomial composition is the reason for PolyXY stuff : */

GtransfoPoly  GtransfoPoly::operator*(const GtransfoPoly &Right) const
{
  // split each transfo into 2d polynomials
  PolyXY plx(*this, 0);
  PolyXY ply(*this, 1);
  PolyXY prx(Right,0);
  PolyXY pry(Right,1);

  // compute the compositions
  PolyXY rx(Composition(plx,prx,pry));
  PolyXY ry(Composition(ply,prx,pry));

  //copy the results the hard way.
  GtransfoPoly result(deg*Right.deg);
  for (unsigned px=0; px<=result.deg; ++px)
    for (unsigned py=0; py<=result.deg-px; ++py)
      {
	result.Coeff(px,py,0) = rx.Coeff(px,py);
	result.Coeff(px,py,1) = ry.Coeff(px,py);
      }
  return result;
}


GtransfoPoly GtransfoPoly::operator+(const GtransfoPoly &Right) const
{
  if (deg >= Right.deg)
    {
      GtransfoPoly res(*this);
      for (unsigned i=0; i<=Right.deg; ++i)
	for (unsigned j = 0; j<=Right.deg-i; ++j)
	  {
	    res.Coeff(i,j,0) += Right.Coeff(i,j,0);
	    res.Coeff(i,j,1) += Right.Coeff(i,j,1);
	  }
      return res;
    }
  else return (Right+(*this));
}

GtransfoPoly GtransfoPoly::operator-(const GtransfoPoly &Right) const
{
  GtransfoPoly res(std::max(deg,Right.deg));
  for (unsigned i=0; i<=res.deg; ++i)
    for (unsigned j = 0; j<=res.deg-i; ++j)
      {
	res.Coeff(i,j,0) = CoeffOrZero(i,j,0) - Right.CoeffOrZero(i,j,0);
	res.Coeff(i,j,1) = CoeffOrZero(i,j,1) - Right.CoeffOrZero(i,j,1);
      }
  return res;
}


void GtransfoPoly::Write(ostream &s) const
{
  s << " GtransfoPoly 1"<< endl;
  s << "degree " << deg << endl;
  int oldprec=s.precision();
  s << setprecision(12);
  for (unsigned k=0;k<2*nterms; ++k)
    s << coeffs[k] << ' ';
  s << endl;
  s << setprecision(oldprec);
}

void GtransfoPoly::Read(istream &s)
{
  int format;
  s >> format;
  if (format != 1)
  throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, " GtransfoPoly::Read : format is not 1 " );

  string degree;
  s >> degree >> deg;
  if (degree != "degree")
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, " GtransfoPoly::Read : expecting \"degree\" and found "+degree );
  SetDegree(deg);
  for (unsigned k=0;k<2*nterms; ++k)
    s >> coeffs[k];
}


std::unique_ptr<GtransfoPoly> InversePolyTransfo(const Gtransfo &Direct, const Frame &F, const double Prec)
{
  StarMatchList sm;
  unsigned nx = 50;
  double stepx = F.Width()/(nx+1);
  unsigned ny = 50;
  double stepy= F.Height()/(ny+1);
  for (unsigned i=0 ;i<nx; ++i)
    for (unsigned j=0; j<ny; ++j)
      {
	Point in((i+0.5)*stepx, (j+0.5)*stepy);
	Point out(Direct.apply(in));
	sm.push_back(StarMatch(out, in , nullptr, nullptr));
      }
  unsigned npairs = sm.size();
  int maxdeg = 9;
  int degree;
  std::unique_ptr<GtransfoPoly> poly;
  for (degree=1; degree<=maxdeg; ++degree)
    {
      poly.reset(new GtransfoPoly(degree));
      poly->fit(sm);
      // compute the chi2 ignoring errors:
      double chi2 = 0;
      for (auto const &i: sm)
        chi2 +=  i.point2.Dist2(poly->apply((i.point1)));
      if (chi2/npairs< Prec*Prec) break;
    }
  if (degree>maxdeg)
    LOGLS_WARN(_log, "InversePolyTransfo: Reached max degree without reaching requested precision: " << Prec);
  return poly;
}



/**************** GtransfoLin ***************************************/
/* GtransfoLin is a specialized constructor of GtransfoPoly
   May be it could just disappear ??
*/

GtransfoLin::GtransfoLin(const double Dx, const double Dy ,
			 const double A11, const double A12,
			 const double A21, const double A22) : GtransfoPoly(1)
{
  dx()  = Dx;
  a11() = A11;
  a12() = A12;
  dy()  = Dy;
  a21() = A21;
  a22() = A22;
}

GtransfoLin::GtransfoLin(const GtransfoPoly &P) : GtransfoPoly(1)
{
  if (P.Degree() !=  1)
      throw pexExcept::InvalidParameterError("Trying to build a GtransfoLin from a higher order transfo. Aborting. ");
  (GtransfoPoly &) (*this) = P;
}


GtransfoLin GtransfoLin::operator*(const  GtransfoLin &Right) const
{
  // There is a general routine in GtransfoPoly that would do the job:
  //  return GtransfoLin(GtransfoPoly::operator*(Right));
  // however, we are using this composition of linear stuff heavily in
  // Gtransfo::LinearApproximation, itself used in InverseTransfo::apply.
  // So, for sake of efficiency, and since it is easy, we take a shortcut:
  GtransfoLin result;
  apply(Right.Dx(), Right.Dy(), result.dx(), result.dy());
  result.a11() = this->A11()*Right.A11() + this->A12()*Right.A21();
  result.a12() = this->A11()*Right.A12() + this->A12()*Right.A22();
  result.a21() = this->A21()*Right.A11() + this->A22()*Right.A21();
  result.a22() = this->A21()*Right.A12() + this->A22()*Right.A22();
  return result;
}

void GtransfoLin::Derivative(const Point &, GtransfoLin &Derivative, const double) const
{
  Derivative = *this;
  Derivative.Coeff(0,0,0) = 0;
  Derivative.Coeff(0,0,1) = 0;
}


GtransfoLin GtransfoLin::LinearApproximation(const Point &, const double) const
{
  return *this;
}

GtransfoLin GtransfoLin::invert() const
{
  //
  //   (T1,M1) * (T2,M2) = 1 i.e (0,1) implies
  //   T1 = -M1 * T2
  //   M1 = M2^(-1)
  //

  double a11 = A11();
  double a12 = A12();
  double a21 = A21();
  double a22 = A22();
  double d = (a11*a22 - a12*a21);
  if (d == 0)
    {
      LOGL_FATAL(_log, "GtransfoLin::invert singular transformation: transfo contents will be dumped to stderr.");
      dump(cerr);
    }

  GtransfoLin result(0,0,a22/d,-a12/d,-a21/d,a11/d);
  double rdx,rdy;
  result.apply(Dx(),Dy(),rdx,rdy);
  result.dx() = -rdx;
  result.dy() = -rdy;
  return result;
}

std::unique_ptr<Gtransfo> GtransfoLin::InverseTransfo(const double, const Frame &) const
{
  return std::unique_ptr<Gtransfo>(new GtransfoLin(this->invert()));
}

double  GtransfoLinRot::fit(const StarMatchList &)
{
  throw pexExcept::NotFoundError("GTransfoLinRot::fit not implemented! aborting");
}


double  GtransfoLinShift::fit(const StarMatchList &List)
{
  int npairs = List.size();
  if (npairs < 3)
    {
      LOGLS_FATAL(_log, "GtransfoLinShift::fit trying to fit a linear transfo with only "
                  << npairs << " matches.");
      return -1;
    }

  double sumr2 = 0; /* used to compute chi2 without relooping */
  /* loop on pairs  and fill */
  Eigen::VectorXd B(2);   B.setZero();
  Eigen::MatrixXd A(2,2);   A.setZero();

  for (auto const &it: List)
    {
      const FatPoint &point1 = it.point1;
      const FatPoint &point2 = it.point2;
      double deltax = point2.x - point1.x;
      double deltay = point2.y - point1.y;
      double vxx = point1.vx+point2.vx;
      double vyy = point1.vy+point2.vy;
      double vxy = point1.vxy+point2.vxy;
      double det = vxx*vyy-vxy*vxy;
      double wxx = vyy/det;
      double wyy = vxx/det;
      double wxy = -vxy/det;
      B(0) +=  deltax*wxx + wxy*deltay;
      B(1) +=  deltay*wyy + wxy*deltax;
      A(0,0) += wxx;
      A(1,1) += wyy;
      A(0,1) += wxy;
      sumr2 += deltax*deltax*wxx + deltay*deltay*wyy
	+2.*wxy*deltax*deltay;
    }
  double det= A(0,0)*A(1,1)-A(0,1)*A(1,0);
  if (det<=0) return -1;
  double tmp = A(0,0);
  A(0,0) = A(1,1)/det;
  A(1,1) = tmp/det;
  A(0,1) = A(1,0) = -A(0,1)/det;
  Eigen::VectorXd sol = A*B;
  (*this) = GtransfoLinShift(sol(0), sol(1));
  return (sumr2  - sol.dot(B)); // chi2 compact form
}


GtransfoLinRot::GtransfoLinRot(const double AngleRad, const Point *Center,
			       const double ScaleFactor)
{
  double c = ScaleFactor*cos(AngleRad);
  double s = ScaleFactor*sin(AngleRad);
  a11() = a22() = c;
  a21() = s;
  a12() = -s;

  // we want that the center does not move : T+M*C = C ==> T = C - M*C
  Point a_point(0.,0.);
  if (Center) a_point = *Center;

  dx() = dy() = 0;
  GtransfoPoly::apply(a_point.x, a_point.y, dx(), dy()); // compute M*C
  dx() = a_point.x - Dx(); dy() = a_point.y - dy();
}


static double deg2rad(double deg)
{
  return deg*M_PI/180.;
}


static double rad2deg(double rad)
{
  return rad*180./M_PI;
}

/*************  WCS transfo ******************/
/************** LinPix2Tan *******************/

/* Implementation note : it seemed wise to incorporate
   the radians to degrees convertion into the linPix2Tan
   part (right in the constructor), and to do the
   opposite operation in the LinPart routine.
   When I was coding the fit, I realized that it was a
   bad idea. then I realized that the fitting routine
   itself was probably useless. Finaly, for a persistor,
   it seems bizarre that the stored data is different
   from what convention (angles in degrees for WCS's)
   would expect.
   So, no more "automatic" degrees to radians and
   radians to degrees conversion. They are explicitely
   done in apply (for TanPix2RaDec and TanRaDec2Pix).
   This is a minor concern though....
*/
BaseTanWcs::BaseTanWcs(const GtransfoLin &Pix2Tan, const Point &TangentPoint,
			const GtransfoPoly *Corrections)
{
  /* the angles returned by linPix2Tan should be in
     degrees. */
  linPix2Tan = Pix2Tan;
  ra0  = deg2rad(TangentPoint.x);
  dec0 = deg2rad(TangentPoint.y);
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = nullptr;
  if (Corrections) corr.reset(new GtransfoPoly(*Corrections));
}

/* with some sort of smart pointer ro handle "corr", we could remove the
   copy constructor, the operator = and the destructor */

// ": Gtransfo" suppresses a warning
BaseTanWcs::BaseTanWcs(const BaseTanWcs &Original) : Gtransfo()
{
  corr = nullptr;
  *this = Original;
}

void  BaseTanWcs::operator = (const BaseTanWcs &Original)
{
  linPix2Tan = Original.linPix2Tan;
  ra0 = Original.ra0;
  dec0 = Original.dec0;
  cos0 = cos(dec0);
  sin0 = sin(dec0);
  corr = nullptr;
  if (Original.corr) corr.reset(new GtransfoPoly(*Original.corr));
}


void BaseTanWcs::apply(const double Xin, const double Yin,
			 double &Xout, double &Yout) const
{
  double l,m; // radians in the tangent plane
  Pix2TP(Xin,Yin,l,m); // l, m in degrees.
  l = deg2rad(l);
  m = deg2rad(m); // now in radians
  // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
  /* At variance with wcslib, it collapses the projection to a plane
     and expression of sidereal cooordinates into a single set of
     operations. */
  double dect = cos0 - m * sin0;
  if (dect == 0)
    {
      LOGL_WARN(_log, "No sidereal coordinates at pole!");
      Xout = 0;
      Yout = 0;
      return;
    }
  double rat = ra0 + atan2(l, dect);
  dect = atan(cos(rat-ra0) * (m * cos0 + sin0) / dect);
  if (rat - ra0 >  M_PI) rat -= (2.*M_PI);
  if (rat - ra0 < -M_PI) rat += (2.*M_PI);
  if (rat < 0.0) rat += (2.*M_PI);
  // convert to deg
  Xout = rad2deg(rat);
  Yout = rad2deg(dect);
}

Point BaseTanWcs::TangentPoint() const
{
  return Point(rad2deg(ra0),rad2deg(dec0));
}

GtransfoLin BaseTanWcs::LinPart() const
{
  return linPix2Tan;
}

void BaseTanWcs::SetCorrections(std::unique_ptr<GtransfoPoly> corrections)
{
    corr = std::move(corrections);
}

Point BaseTanWcs::CrPix() const
{
  /* CRPIX's are defined by:
                    ( CD1_1  CD1_2 )   (X - crpix1)
     transformed =  (              ) * (          )
                    ( CD2_1  CD2_2 )   (Y - crpix2)

     so that CrPix is the point which transforms to (0,0)
  */
  const GtransfoLin inverse = linPix2Tan.invert();
  return Point(inverse.Dx(),inverse.Dy());
}



BaseTanWcs::~BaseTanWcs() { }

/*************************** TanPix2RaDec ***************/

TanPix2RaDec::TanPix2RaDec(const GtransfoLin &Pix2Tan,
			   const Point &TangentPoint,
			   const GtransfoPoly* Corrections) : BaseTanWcs(Pix2Tan, TangentPoint, Corrections)
{
}


// ": Gtransfo" suppresses a warning
TanPix2RaDec::TanPix2RaDec() : BaseTanWcs(GtransfoLin(), Point(0,0), nullptr)
{
}


std::unique_ptr<Gtransfo> TanPix2RaDec::ReduceCompo(const Gtransfo *Right) const
{
  const GtransfoLin *lin = dynamic_cast<const GtransfoLin *>(Right);
  if (lin && lin->Degree() == 1) return std::unique_ptr<Gtransfo>(new TanPix2RaDec((*this)*(*lin)));
  return std::unique_ptr<Gtransfo>(nullptr);
}


TanPix2RaDec TanPix2RaDec::operator *(const GtransfoLin &Right) const
{
  TanPix2RaDec result(*this);
  result.linPix2Tan = result.linPix2Tan * Right;
  return result;
}

TanRaDec2Pix TanPix2RaDec::invert() const
{
  if (corr != nullptr)
    {
      LOGL_WARN(_log, "You are inverting a TanPix2RaDec with corrections.");
      LOGL_WARN(_log, "The inverse you get ignores the corrections!");
    }
  return TanRaDec2Pix(LinPart().invert(),TangentPoint());
}

std::unique_ptr<Gtransfo> TanPix2RaDec::RoughInverse(const Frame &) const
{
  return std::unique_ptr<Gtransfo>(new TanRaDec2Pix(LinPart().invert(),TangentPoint()));
}

std::unique_ptr<Gtransfo> TanPix2RaDec::InverseTransfo(const double Precision,
					const Frame& Region) const
{
  if (!corr) return std::unique_ptr<Gtransfo>(new TanRaDec2Pix(LinPart().invert(),TangentPoint()));
  else return std::unique_ptr<Gtransfo>(new GtransfoInverse(this, Precision, Region));
}


GtransfoPoly TanPix2RaDec::Pix2TangentPlane() const
{
  if (corr) return (*corr)*linPix2Tan;
  else return linPix2Tan;
}

void TanPix2RaDec::Pix2TP(double Xin, double Yin,
			  double &Xtp, double & Ytp) const
{
  linPix2Tan.apply(Xin, Yin, Xtp, Ytp); // Xtp, Ytp in degrees.
  if (corr)
    {
      double xtmp = Xtp;
      double ytmp = Ytp;
      corr->apply(xtmp , ytmp, Xtp, Ytp); // still in degrees.
    }
}


std::unique_ptr<Gtransfo> TanPix2RaDec::Clone() const
{
  return std::unique_ptr<Gtransfo>(new TanPix2RaDec(LinPart(), TangentPoint(), corr.get()));
}

void TanPix2RaDec::dump(ostream &stream) const
{
  stream << " TanPix2RaDec, lin part :" << endl << linPix2Tan;
  Point tp = TangentPoint();
  stream << " tangent point " << tp.x << ' ' << tp.y << endl;
  Point crpix = CrPix();
  stream << " crpix " << crpix.x << ' ' << crpix.y << endl;
  if (corr) stream << "PV correction: " << endl << *corr;
}


double  TanPix2RaDec::fit(const StarMatchList &)
{
  /* OK we could implement this routine, but it is
     probably useless since to do the match, we have to
     project from sky to tangent plane. When a match is
     found, it is easier to carry out the fit in the
     tangent plane, rather than going back to the celestial
     sphere (and reproject to fit...). Anyway if this
     message shows up, we'll think about it.
  */
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "TanPix2RaDec::fit is NOT implemented (although it is doable)) ");
  return -1;
}

/*************************** TanSipPix2RaDec ***************/

TanSipPix2RaDec::TanSipPix2RaDec(const GtransfoLin &Pix2Tan,
			   const Point &TangentPoint,
			   const GtransfoPoly* Corrections) : BaseTanWcs(Pix2Tan, TangentPoint, Corrections)
{
}


// ": Gtransfo" suppresses a warning
TanSipPix2RaDec::TanSipPix2RaDec() : BaseTanWcs(GtransfoLin(), Point(0,0), nullptr)
{
}



/* Would require some checks before cooking up something more efficient
   than just a linear approximation */
#if 0
std::unique_ptr<Gtransfo> TanPix2RaDec::RoughInverse(const Frame &Region) const
{
  if (&Region) {}
  return std::unique_ptr<Gtransfo>(new TanRaDec2Pix(LinPart().invert(),TangentPoint()));
}
#endif

std::unique_ptr<Gtransfo> TanSipPix2RaDec::InverseTransfo(const double Precision,
					const Frame& Region) const
{/* We have not implemented (yet) the reverse corrections available in SIP */
   return std::unique_ptr<Gtransfo>(new GtransfoInverse(this, Precision, Region));
}


GtransfoPoly TanSipPix2RaDec::Pix2TangentPlane() const
{
  if (corr) return GtransfoPoly(linPix2Tan)*(*corr);
  else return linPix2Tan;
}

void TanSipPix2RaDec::Pix2TP(double Xin, double Yin, double &Xtp, double &Ytp) const
{ // Xtp, Ytp returned in degrees
  if (corr)
    {
      double xtmp, ytmp;
      corr->apply(Xin, Yin, xtmp, ytmp);
      linPix2Tan.apply(xtmp, ytmp, Xtp, Ytp);
    }
  else linPix2Tan.apply(Xin, Yin, Xtp, Ytp);
}


std::unique_ptr<Gtransfo> TanSipPix2RaDec::Clone() const
{
  return std::unique_ptr<Gtransfo>(new TanSipPix2RaDec(LinPart(), TangentPoint(), corr.get()));
}

void TanSipPix2RaDec::dump(ostream &stream) const
{
  stream << " TanSipPix2RaDec, lin part :" << endl << linPix2Tan;
  Point tp = TangentPoint();
  stream << " tangent point " << tp.x << ' ' << tp.y << endl;
  Point crpix = CrPix();
  stream << " crpix " << crpix.x << ' ' << crpix.y << endl;
  if (corr) stream << "PV correction: " << endl << *corr;
}


double  TanSipPix2RaDec::fit(const StarMatchList &)
{
  /* OK we could implement this routine, but it is
     probably useless since to do the match, we have to
     project from sky to tangent plane. When a match is
     found, it is easier to carry out the fit in the
     tangent plane, rather than going back to the celestial
     sphere (and reproject to fit...). Anyway if this
     message shows up, we'll think about it.
  */
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "TanSipPix2RaDec::fit is NOT implemented (although it is doable)) ");
  return -1;
}



/***************  reverse transfo of TanPix2RaDec: TanRaDec2Pix ********/


TanRaDec2Pix::TanRaDec2Pix(const GtransfoLin &Tan2Pix, const Point &TangentPoint) : linTan2Pix(Tan2Pix)
{
  SetTangentPoint(TangentPoint);
}

void TanRaDec2Pix::SetTangentPoint(const Point &TangentPoint)
{
/* the radian to degrees conversion after projection
    is handled in apply */
  ra0  = deg2rad(TangentPoint.x);
  dec0 = deg2rad(TangentPoint.y);
  cos0 = cos(dec0);
  sin0 = sin(dec0);
}


TanRaDec2Pix::TanRaDec2Pix() : linTan2Pix()
{
  ra0 = dec0 = 0;
  cos0 = 1;
  sin0 = 0;
}



Point TanRaDec2Pix::TangentPoint() const
{
  return Point(rad2deg(ra0),rad2deg(dec0));
}

GtransfoLin TanRaDec2Pix::LinPart() const
{
  return linTan2Pix;
}

// Use analytic derivatives, computed at the same time as the transform itself
void TanRaDec2Pix::TransformPosAndErrors(const FatPoint &In,
					 FatPoint &Out) const
{
  /* this routine is very similar to apply, but also propagates errors.
     The deg2rad and rad2deg are ignored for errors because they act as
     2 global scalings that cancel each other.
     Derivatives were computed using maple:

     l1 := sin(a - a0)*cos(d);
     m1 := sin(d)*sin(d0)+cos(d)*cos(d0)*cos(a-a0);
     l2 := sin(d)*cos(d0)-cos(d)*sin(d0)*cos(a-a0);
     simplify(diff(l1/m1,a));
     simplify(diff(l1/m1,d));
     simplify(diff(l2/m1,a));
     simplify(diff(l2/m1,d));

     Checked against Gtransfo::TransformPosAndErrors (dec 09)
  */
  double ra = deg2rad(In.x);
  double dec = deg2rad(In.y);
  if (ra-ra0 >  M_PI ) ra -= (2.* M_PI);
  if (ra-ra0 < -M_PI ) ra += (2.* M_PI);
  // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
  // The same code is copied in ::apply()

  double coss = cos(dec);
  double sins = sin(dec);
  double sinda = sin(ra-ra0);
  double cosda = cos(ra -ra0);
  double l = sinda * coss;
  double m = sins * sin0 + coss * cos0 * cosda;
  l = l/m;
  m = (sins * cos0 - coss * sin0 * cosda)/ m;

  // derivatives
  double deno = sq(sin0)-sq(coss)+sq(coss*cos0)*(1+sq(cosda))+2*sins*sin0*coss*cos0*cosda;
  double a11 = coss*(cosda*sins*sin0+coss*cos0)/deno;
  double a12 = -sinda*sin0/deno;
  double a21 = coss*sinda*sins/deno;
  double a22 = cosda/deno;

  FatPoint tmp;
  tmp.vx = a11*(a11*In.vx + 2*a12*In.vxy) + a12*a12*In.vy;
  tmp.vy = a21*a21*In.vx + a22*a22*In.vy + 2.*a21*a22*In.vxy;
  tmp.vxy = a21*a11*In.vx + a22*a12*In.vy + (a21*a12+a11*a22)*In.vxy;

  // l and m are now coordinates in the tangent plane, in radians.
  tmp.x = rad2deg(l);
  tmp.y = rad2deg(m);

  linTan2Pix.TransformPosAndErrors(tmp, Out);
}


void TanRaDec2Pix::apply(const double Xin, const double Yin, double &Xout, double &Yout) const
{
  double ra = deg2rad(Xin);
  double dec = deg2rad(Yin);
  if (ra-ra0 >  M_PI ) ra -= (2.* M_PI);
  if (ra-ra0 < -M_PI ) ra += (2.* M_PI);
  // Code inspired from worldpos.c in wcssubs (ancestor of the wcslib)
  // The same code is copied in ::TransformPosAndErrors()
  double coss = cos(dec);
  double sins = sin(dec);
  double l = sin(ra-ra0) * coss;
  double m = sins * sin0 + coss * cos0 * cos (ra -ra0);
  l = l/m;
  m = (sins * cos0 - coss * sin0 * cos( ra -ra0))/ m;
  // l and m are now coordinates in the tangent plane, in radians.
  l = rad2deg(l);
  m = rad2deg(m);
  linTan2Pix. apply(l,m, Xout, Yout);
}


TanPix2RaDec TanRaDec2Pix::invert() const
{
  return TanPix2RaDec(LinPart().invert(),TangentPoint());
}


void TanRaDec2Pix::dump(ostream &stream) const
{
  Point tp = TangentPoint();
  stream << " tan2pix " << linTan2Pix << " tangent point " << tp.x << ' ' << tp.y << endl;
}

std::unique_ptr<Gtransfo> TanRaDec2Pix::RoughInverse(const Frame &) const
{
  return std::unique_ptr<Gtransfo>(new TanPix2RaDec(LinPart().invert(),TangentPoint()));
}

std::unique_ptr<Gtransfo> TanRaDec2Pix::InverseTransfo(const double, const Frame &) const
{
  return std::unique_ptr<Gtransfo>(new TanPix2RaDec(LinPart().invert(),TangentPoint()));
}


std::unique_ptr<Gtransfo> TanRaDec2Pix::Clone() const
{
  return std::unique_ptr<Gtransfo>(new TanRaDec2Pix(*this));
}

double TanRaDec2Pix::fit(const StarMatchList &)
{
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "TanRaDec2Pix::fit is NOT implemented (although it is doable)) ");
  return -1;
}

/*************  a "run-time" transfo, that does not require to
modify this file */

UserTransfo::UserTransfo(GtransfoFun &Fun, const void *UserData):
  userFun(Fun), userData(UserData)
{}

void UserTransfo::apply(const double Xin, const double Yin,
			double &Xout, double &Yout) const
{
  userFun(Xin,Yin,Xout,Yout, userData);
}

void UserTransfo::dump(ostream &stream) const
{
  stream << "UserTransfo with user function @ " << userFun
	 << "and userData@ " << userData << endl;
}

double UserTransfo::fit(const StarMatchList &)
{
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "UserTransfo::fit is NOT implemented (and will never be)) ");
  return -1;
}

std::unique_ptr<Gtransfo> UserTransfo::Clone() const
{
  return std::unique_ptr<Gtransfo>(new UserTransfo(*this));
}


/*************************************************************/


std::unique_ptr<Gtransfo> GtransfoRead(const std::string &FileName)
{
  ifstream s(FileName.c_str());
  if (!s)
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, " GtransfoRead : cannot open " + FileName);
  try
    {
        std::unique_ptr<Gtransfo> res(GtransfoRead(s));
      s.close();
      return res;
    }
  catch (pex::exceptions::InvalidParameterError e)
    {
      throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
			std::string(e.what())+" in file "+FileName);
    }
}

std::unique_ptr<Gtransfo> GtransfoRead(istream &s)
{
  std::string type;
  s >> type;
  if (s.fail())
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "GtransfoRead : could not find a Gtransfotype");
  if (type == "GtransfoIdentity")
    {std::unique_ptr<GtransfoIdentity> res(new GtransfoIdentity()); res->Read(s); return std::move(res);}
  else if (type == "GtransfoPoly")
    {std::unique_ptr<GtransfoPoly> res(new GtransfoPoly()); res->Read(s); return std::move(res);}
  else
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, " GtransfoRead : No reader for Gtransfo type "+ type);
}


}} // end of namespaces
