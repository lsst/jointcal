#include <iostream>
#include <iomanip>
#include <algorithm>
#include "lsst/jointcal/AstromFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Mapping.h"

#include "lsst/jointcal/Gtransfo.h"
#include "Eigen/Sparse"
#include "Eigen/CholmodSupport" // to switch to cholmod
#include <time.h> // for clock
#include "lsst/pex/exceptions.h"
#include <fstream>
#include "lsst/jointcal/Tripletlist.h"

typedef Eigen::SparseMatrix<double> SpMat;


#if 0
template<typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSupernodalLLT2 : public Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSupernodalLLT2<_MatrixType, _UpLo> >
{
  typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSupernodalLLT2> Base;
    using Base::m_cholmod;
    
  public:
    
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;
    
    CholmodSupernodalLLT2() : Base() { init(); }
  //! The factorization happens in the constructor.
    CholmodSupernodalLLT2(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

  //! this routine is the one we added.
  int update(const SpMat &H, const bool UpOrDown)
  {
    // check size
    const Index size = Base::m_cholmodFactor->n;
    EIGEN_UNUSED_VARIABLE(size);
    eigen_assert(size==H.rows());
    cholmod_sparse C_cs = viewAsCholmod(H);
    /* We have to apply the magic permutation to the update matrix, 
       read page 117 of Cholmod UserGuide.pdf */
    cholmod_sparse *C_cs_perm = cholmod_submatrix(C_cs,
						  Base::m_cholmodFactor->Perm,
						  Base::m_cholmodFactor->n,
						  NULL, -1, NULL, true, true,
						  &this->cholmod());
    
    int ret = cholmod_updown(UpOrDown, &C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
    cholmod_free_sparse(C_cs_perm,  &this->cholmod());
    return ret;
  }

    ~CholmodSupernodalLLT2() {}
  protected:
    void init()
    {
      m_cholmod.final_asis = 1;
      m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
      // In CholmodBase::CholmodBase(), the following statement is missing in
      // SuiteSparse 3.2.0.8. Fixed in 3.2.7
      Base::m_shiftOffset[0] = Base::m_shiftOffset[1] = RealScalar(0.0);
    }
};
#endif


//! Cholesky factorization class using cholmod, with the small-rank update capability. 
/*! Class derived from Eigen's CholmodBase, to add the factorization
    update capability to the interface. Besides this addition, it
    behaves the same way as Eigen's native Cholesky factorization
    classes. It relies on the simplicial LDLt factorization.*/
template<typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSimplicialLDLT2 : public Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2<_MatrixType, _UpLo> >
{
  typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2> Base;
    using Base::m_cholmod;
    
  public:
    
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;
    typedef typename MatrixType::RealScalar RealScalar;
    
    CholmodSimplicialLDLT2() : Base() { init(); }

    CholmodSimplicialLDLT2(const MatrixType& matrix) : Base()
    {
      init();
      this->compute(matrix);
    }

    // this routine is the one we added
    int update(const SpMat &H, const bool UpOrDown)
    {
      // check size
      const Index size = Base::m_cholmodFactor->n;
      EIGEN_UNUSED_VARIABLE(size);
      eigen_assert(size==H.rows());
      
      cholmod_sparse C_cs = viewAsCholmod(H);
      /* We have to apply the magic permutation to the update matrix, 
	 read page 117 of Cholmod UserGuide.pdf */
      cholmod_sparse *C_cs_perm = cholmod_submatrix(&C_cs,
						    (int *) Base::m_cholmodFactor->Perm,
						    Base::m_cholmodFactor->n,
						    NULL, -1, true, true,
						    &this->cholmod());
      assert(C_cs_perm);
      int ret = cholmod_updown(UpOrDown, C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
      cholmod_free_sparse(&C_cs_perm,  &this->cholmod());
      assert(ret != 0);
      return ret;
    }

    ~CholmodSimplicialLDLT2() {}
  protected:
    void init()
    {
      m_cholmod.final_asis = 1;
      m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
      // In CholmodBase::CholmodBase(), the following statement is missing in
      // SuiteSparse 3.2.0.8. Fixed in 3.2.7
      Base::m_shiftOffset[0] = Base::m_shiftOffset[1] = RealScalar(0.0);
    }
};



using namespace std;

static double sqr(const double &x) {return x*x;}

//const double posErrorIncrement=0.02;

namespace lsst {
namespace jointcal {


AstromFit::AstromFit(Associations &A, DistortionModel *D, double PosError) : 
  _assoc(A),  _distortionModel(D), _posError(PosError)
{
  _LastNTrip = 0;
  _JDRef = 0;

  _posError = PosError;

  _referenceColor = 0; 
  _sigCol = 0;
  unsigned count = 0;
  for (auto i=_assoc.fittedStarList.begin(); 
       i!= _assoc.fittedStarList.end() ; ++i) 
    {
      _referenceColor += (*i)->color;
      _sigCol += sqr((*i)->color);
      count++;
    }
  if (count)
    {
      _referenceColor /= double(count);
      if (_sigCol>0) _sigCol = sqrt(_sigCol/count - sqr(_referenceColor));
    }
  cout << "INFO: reference Color : " << _referenceColor << " sig " << _sigCol << endl;

  _nRefrac = _assoc.NBands();
  _refracCoefficient.resize(_nRefrac,0);

  _nMeasuredStars = 0;
  // The various _npar... are initialized in AssignIndices.
  // Although there is no reason to adress them before one might be tempted by
  // evaluating a Chi2 rightaway, .. which uses these counts, so:
  AssignIndices("");

}



#define NPAR_PM 2

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstromFit::TransformFittedStar(const FittedStar &F,
				     const Gtransfo * Sky2TP,
				     const Point &RefractionVector,
				     const double &RefractionCoeff,
				     const double &Jd) const
{
  Point fittedStarInTP =  Sky2TP->apply(F);
  if (F.mightMove)
    {
      fittedStarInTP.x += F.pmx*Jd;
      fittedStarInTP.y += F.pmy*Jd;
    }
  // account for atmospheric refraction: does nothing if color 
  // have not been assigned
  // the color definition shouldbe the same when computing derivatives
  double color = F.color - _referenceColor;
  fittedStarInTP.x += RefractionVector.x * color * RefractionCoeff;
  fittedStarInTP.y += RefractionVector.y * color * RefractionCoeff;
  return fittedStarInTP;
}

/*! This is the first implementation of an error "model".  We'll
  certainly have to upgrade it. MeasuredStar provides the mag in case
  we need it.  */
static void TweakAstromMeasurementErrors(FatPoint &P, const MeasuredStar &Ms, double error)
{
  static bool called=false;
  static double increment = 0;
  if (!called)
    {
      increment = sqr(error); // was in Preferences
      called = true;
    }
  P.vx += increment;
  P.vy += increment;
}

static bool heavyDebug = false;
static unsigned fsIndexDebug = 0;



// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstromFit::LSDerivatives1(const CcdImage &Ccd, 
			      TripletList &TList, Eigen::VectorXd &Rhs,
			      const MeasuredStarList *M) const
{
  /***************************************************************************/
  /**  Changes in this routine should be reflected into AccumulateStatImage  */
  /***************************************************************************/
  /* Setup */
  /* this routine works in two different ways: either providing the
     Ccd, of providing the MeasuredStarList. In the latter case, the
     Ccd should match the one(s) in the list. */
  if (M) 
    assert ( (*(M->begin()))->ccdImage == &Ccd);

  // get the Mapping
  const Mapping *mapping = _distortionModel->GetMapping(Ccd);
  // count parameters
  unsigned npar_mapping = (_fittingDistortions) ? mapping->Npar() : 0;
  unsigned npar_pos = (_fittingPos) ? 2 : 0;
  unsigned npar_refrac = (_fittingRefrac) ? 1 : 0;
  unsigned npar_pm = (_fittingPM) ? NPAR_PM : 0;
  unsigned npar_tot =  npar_mapping + npar_pos + npar_refrac + npar_pm;
  // if (npar_tot == 0) this CcdImage does not contribute 
  // any constraint to the fit, so :
  if (npar_tot == 0) return;
  vector<unsigned> indices(npar_tot,-1);
  if (_fittingDistortions)  mapping->GetMappingIndices(indices);

  // proper motion stuff
  double jd = Ccd.JD() - _JDRef;
  // refraction stuff
  Point refractionVector = Ccd.ParallacticVector();
  unsigned iband = Ccd.BandRank();
  double refractionCoefficient = _refracCoefficient.at(iband);
  // transformation from sky to TP
  const Gtransfo* sky2TP = _distortionModel->Sky2TP(Ccd);
  // reserve matrices once for all measurements
  GtransfoLin dypdy;
  // the shape of h (et al) is required this way in order to be able to 
  // separate derivatives along x and y as vectors.
  Eigen::MatrixX2d h(npar_tot,2), halpha(npar_tot,2), hw(npar_tot,2); 
  Eigen::Matrix2d transW(2,2);
  Eigen::Matrix2d alpha(2,2);
  Eigen::VectorXd grad(npar_tot);
  // current position in the Jacobian
  unsigned kTriplets = TList.NextFreeIndex();
  const MeasuredStarList &catalog = (M) ? *M : Ccd.CatalogForFit();

  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      const MeasuredStar& ms = **i;
      if (!ms.IsValid()) continue;
      // tweak the measurement errors
      FatPoint inPos = ms;
      TweakAstromMeasurementErrors(inPos, ms, _posError);
      h.setZero(); // we cannot be sure that all entries will be overwritten.
      FatPoint outPos;
      // should *not* fill h if WhatToFit excludes mapping parameters.
      if (_fittingDistortions) 
	  mapping->ComputeTransformAndDerivatives(inPos, outPos, h);
      else mapping->TransformPosAndErrors(inPos,outPos);

      unsigned ipar = npar_mapping;
      double det = outPos.vx*outPos.vy-sqr(outPos.vxy);
      if (det <=0 || outPos.vx <=0 || outPos.vy<=0) {
	cout << " WARNING: inconsistent measurement errors :drop measurement at " << Point(ms) << " in image " << Ccd.Name() << endl;
	continue;
      }	
      transW(0,0) = outPos.vy/det;
      transW(1,1) = outPos.vx/det;
      transW(0,1) = transW(1,0) = -outPos.vxy/det;
      // compute alpha, a triangular square root
      // of transW (i.e. a Cholesky factor)
      alpha(0,0) = sqrt(transW(0,0));
      // checked that  alpha*alphaT = transW
      alpha(1,0) = transW(0,1)/alpha(0,0); 
      // DB - I think that the next line is equivalent to : alpha(1,1) = 1./sqrt(outPos.vy)
      // PA - seems correct !
      alpha(1,1) = 1./sqrt(det*transW(0,0));
      alpha(0,1) = 0;
      
      const FittedStar *fs = ms.GetFittedStar();

      Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						 refractionVector, 
						 refractionCoefficient,
						 jd);

      // compute derivative of TP position w.r.t sky position ....
      if (npar_pos>0) // ... if actually fitting FittedStar position
	{
	  sky2TP->Derivative(*fs, dypdy, 1e-3);
	  // sign checked
	  // TODO Still have to check with non trivial non-diagonal terms
	  h(npar_mapping,0) = -dypdy.A11();
	  h(npar_mapping+1,0) = -dypdy.A12();
	  h(npar_mapping,1) = -dypdy.A21();	  
	  h(npar_mapping+1,1) = -dypdy.A22();
	  indices[npar_mapping] = fs->IndexInMatrix();
	  indices.at(npar_mapping+1) = fs->IndexInMatrix()+1;
	  ipar += npar_pos;
	}
      /* only consider proper motions of objects allowed to move,
	 unless the fit is going to be degenerate */
      if (_fittingPM && fs->mightMove)
	{
	  h(ipar,0) = -jd; // Sign unchecked but consistent with above
	  h(ipar+1, 1) = -jd;
	  indices[ipar] = fs->IndexInMatrix()+2;
	  indices[ipar+1] = fs->IndexInMatrix()+3;
	  ipar+= npar_pm;
	}
      if (_fittingRefrac)
	{
	  /* if the definition of color changes, it has to remain
	     consistent with TransformFittedStar */
	  double color = fs->color - _referenceColor;
	  // sign checked 
	  h(ipar,0) = -refractionVector.x*color;
	  h(ipar,1) = -refractionVector.y*color;
	  indices[ipar] = _refracPosInMatrix+iband;
	  ipar += 1;
	}

      // We can now compute the residual
      Eigen::Vector2d res(fittedStarInTP.x-outPos.x, fittedStarInTP.y-outPos.y);

      // do not write grad = h*transW*res to avoid 
      // dynamic allocation of a temporary
      halpha = h*alpha;
      hw = h*transW;
      grad = hw*res;
      // now feed in triplets and Rhs
      for (unsigned ipar=0; ipar<npar_tot; ++ipar)
	{
	  for (unsigned  ic=0; ic<2; ++ic)
	    {
	      double val = halpha(ipar,ic);
	      if (val ==0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
	      TList.AddTriplet(indices[ipar], kTriplets+ic,val);
#else
	      TList.AddTriplet(kTriplets+ic, indices[ipar], val);
#endif
	    }
	  Rhs(indices[ipar]) += grad(ipar); 
	}
      kTriplets += 2; // each measurement contributes 2 columns in the Jacobian
    } // end loop on measurements
  TList.SetNextFreeIndex(kTriplets);
}    

// we could consider computing the chi2 here.
// (although it is not extremely useful)

#define HACK_REF_ERRORS 1. // used to isolate the measurement or ref terms

void AstromFit::LSDerivatives2(const FittedStarList &Fsl, TripletList &TList, Eigen::VectorXd &Rhs) const
{
  /* We compute here the derivatives of the terms involving fitted
     stars and reference stars. They only provide contributions if we
     are fitting positions: */
  if (! _fittingPos) return;
  /* the other case where the accumulation of derivatives stops 
     here is when there are no RefStars */
  if (_assoc.refStarList.size() == 0) return;
  Eigen::Matrix2d w(2,2);
  Eigen::Matrix2d alpha(2,2);
  Eigen::Matrix2d h(2,2), halpha(2,2), hw(2,2);
  GtransfoLin der;
  Eigen::Vector2d res,grad;
  unsigned indices[2+NPAR_PM];
  unsigned kTriplets = TList.NextFreeIndex();
  /* We cannot use the spherical coordinates directly to evaluate
     Euclidean distances, we have to use a projector on some plane in
     order to express least squares. Not projecting could lead to a
     disaster around the poles or across alpha=0.  So we need a
     projector. We construct a projector and will change its
     projection point at every object */
  TanRaDec2Pix proj(GtransfoLin(), Point(0.,0.));
  for (auto i = Fsl.cbegin(); i!= Fsl.end(); ++i)
    {
      const FittedStar &fs = **i;
      const RefStar *rs = fs.GetRefStar();
      if (rs == NULL) continue;
      proj.SetTangentPoint(fs);
      // fs projects to (0,0), no need to compute its transform.
      FatPoint rsProj;
      proj.TransformPosAndErrors(*rs, rsProj);
      proj.Derivative(fs, der, 1e-4);
      // sign checked. TODO check that the off-diagonal terms are OK.
      h(0,0) = -der.A11();
      h(1,0) = -der.A12();
      h(0,1) = -der.A21();
      h(1,1) = -der.A22();
      // TO DO : account for proper motions.
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      if (rsProj.vx <=0 || rsProj.vy <=0 || det <= 0) 
	{
	  cout << " WARNING: Ref star error matrix not posdef:  " << endl
	       << *rs << endl;
	  continue;
	}
      w(0,0) = rsProj.vy/det;
      w(0,1) = w(1,0) = -rsProj.vxy/det;
      w(1,1) = rsProj.vx/det;
      w *= HACK_REF_ERRORS;
      det /= sqr(HACK_REF_ERRORS);
      // compute alpha, a triangular square root
      // of w (i.e. a Cholesky factor)
      alpha(0,0) = sqrt(w(0,0));
      // checked that  alpha*alphaT = transW
      alpha(1,0) = w(0,1)/alpha(0,0); 
      alpha(1,1) = 1./sqrt(det*w(0,0));
      alpha(0,1) = 0;
      indices[0] = fs.IndexInMatrix();
      indices[1] = fs.IndexInMatrix()+1;
      unsigned npar_tot = 2;
      /* TODO: account here for proper motions in the reference
      catalog. We can code the effect and set the value to 0. Most
      (all?)  catalogs do not even come with a reference epoch. Gaia
      will change that. When refraction enters into the game, one should 
      pay attention to the orientation of the frame */

      /* The residual should be Proj(fs)-Proj(*rs) in order to be consistent
	 with the measurement terms. Since P(fs) = 0, we have: */
      res[0] = -rsProj.x;
      res[1] = -rsProj.y;
      halpha = h*alpha;
      // grad = h*w*res
      hw = h*w;
      grad = hw*res;
      // now feed in triplets and Rhs
      for (unsigned ipar=0; ipar<npar_tot; ++ipar)
	{
	  for (unsigned  ic=0; ic<2; ++ic)
	    {
	      double val = halpha(ipar,ic);
	      if (val ==0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
	      TList.AddTriplet(indices[ipar], kTriplets+ic,val);
#else
	      TList.AddTriplet(kTriplets+ic, indices[ipar], val);
#endif
	    }
	  Rhs(indices[ipar]) += grad(ipar); 
	}
      kTriplets += 2; // each measurement contributes 2 columns in the Jacobian      
    }
  TList.SetNextFreeIndex(kTriplets);
}


//! this routine computes the derivatives of all LS terms, including the ones that refer to references stars, if any
void AstromFit::LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs) const
{
  auto L = _assoc.TheCcdImageList();
  for (auto im=L.cbegin(); im!=L.end() ; ++im)
    {
      LSDerivatives1(**im, TList, Rhs);
    }
  LSDerivatives2(_assoc.fittedStarList, TList, Rhs);
}


// This is almost a selection of lines of LSDerivatives1(CcdImage ...)
/* This routine (and the following one) is template because it is used
both with its first argument as "const CCdImage &" and "CcdImage &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ImType, class Accum> 
void AstromFit::AccumulateStatImage(ImType &Ccd, Accum &Accu) const
{
  /**********************************************************************/
  /**  Changes in this routine should be reflected into LSDerivatives1  */
  /**********************************************************************/
  /* Setup */
  // 1 : get the Mapping's
  const Mapping *mapping = _distortionModel->GetMapping(Ccd);
  // proper motion stuff
  double jd = Ccd.JD() - _JDRef;
  // refraction stuff
  Point refractionVector = Ccd.ParallacticVector();
  double refractionCoefficient = _refracCoefficient.at(Ccd.BandRank());
  // transformation from sky to TP
  const Gtransfo* sky2TP = _distortionModel->Sky2TP(Ccd);
  // reserve matrix once for all measurements
  Eigen::Matrix2Xd transW(2,2);

  auto &catalog = Ccd.CatalogForFit();
  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      auto &ms = **i;
      if (!ms.IsValid()) continue;
      // tweak the measurement errors
      FatPoint inPos = ms;
      TweakAstromMeasurementErrors(inPos, ms, _posError);

      FatPoint outPos;
      // should *not* fill h if WhatToFit excludes mapping parameters.
      mapping->TransformPosAndErrors(inPos, outPos);
      double det = outPos.vx*outPos.vy-sqr(outPos.vxy);
      if (det <=0 || outPos.vx <=0 || outPos.vy<=0) {
	cout << " WARNING: inconsistent measurement errors :drop measurement at " << Point(ms) << " in image " << Ccd.Name() << endl;
	continue;
      }	
      transW(0,0) = outPos.vy/det;
      transW(1,1) = outPos.vx/det;
      transW(0,1) = transW(1,0) = -outPos.vxy/det;

      const FittedStar *fs = ms.GetFittedStar();
      Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						 refractionVector, 
						 refractionCoefficient,
						 jd);

      Eigen::Vector2d res(fittedStarInTP.x-outPos.x, fittedStarInTP.y-outPos.y); 
      double chi2Val = res.transpose()*transW*res;

      Accu.AddEntry(chi2Val, 2, &ms);
    }// end of loop on measurements
}


//! for a list of images.
template <class ListType, class Accum> 
void AstromFit::AccumulateStatImageList(ListType &L, Accum &Accu) const
{
  for (auto im=L.begin(); im!=L.end() ; ++im)
    {
      AccumulateStatImage(**im, Accu);
    }
}

template <class Accum> 
void AstromFit::AccumulateStatRefStars(Accum &Accu) const
{
  /* If you wonder why we project here, read comments in 
     AstromFit::LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) */
  FittedStarList &fsl = _assoc.fittedStarList;
  TanRaDec2Pix proj(GtransfoLin(), Point(0.,0.));
  for (auto i = fsl.begin(); i!= fsl.end(); ++i)
    {
      FittedStar &fs = **i;
      const RefStar *rs = fs.GetRefStar();
      if (rs == NULL) continue;
      proj.SetTangentPoint(fs);
      // fs projects to (0,0), no need to compute its transform.
      FatPoint rsProj;
      proj.TransformPosAndErrors(*rs, rsProj);
      // TO DO : account for proper motions.
      double rx = rsProj.x; // -fsProj.x (which is 0)
      double ry = rsProj.y;
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      double wxx = rsProj.vy/det;
      double wyy = rsProj.vx/det;
      double wxy = -rsProj.vxy/det;
      wxx *= HACK_REF_ERRORS;
      wyy *= HACK_REF_ERRORS;
      wxy *= HACK_REF_ERRORS;
      Accu.AddEntry(wxx*sqr(rx) + 2*wxy*rx*ry+ wyy*sqr(ry), 2, &fs);
    }
}


//! for the list of images in the provided  association and the reference stars, if any
Chi2 AstromFit::ComputeChi2() const
{
  Chi2 chi2;
  AccumulateStatImageList(_assoc.TheCcdImageList(), chi2);
  // now ref stars:
  AccumulateStatRefStars(chi2);
  // so far, ndof contains the number of squares.
  // So, subtract here the number of parameters.
  chi2.ndof -= _nParTot;
  return chi2;
}

//! a class to accumulate chi2 contributions together with pointers to the contributors.
/*! This structure allows to compute the chi2 statistics (average and
  variance) and directly point back to the bad guys without
  relooping. The Chi2Entry routine makes it compatible with
  AccumulateStatImage and AccumulateStatImageList. */
struct Chi2Entry
{
  double chi2;
  BaseStar *ps;

  Chi2Entry(const double &c, BaseStar *s): chi2(c), ps(s) {}
  // for sort
  bool operator < (const Chi2Entry &R) const {return (chi2<R.chi2);}
};

struct Chi2Vect : public std::vector<Chi2Entry>
{
  void AddEntry(const double &Chi2Val, unsigned ndof, BaseStar *ps)
  { this->push_back(Chi2Entry(Chi2Val,ps));}

};

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void AstromFit::GetMeasuredStarIndices(const MeasuredStar &Ms, 
				       std::vector<unsigned> &Indices) const
{
  if (_fittingDistortions)
    {
      const Mapping *mapping = _distortionModel->GetMapping(*Ms.ccdImage);
      mapping->GetMappingIndices(Indices);
    }
  const FittedStar *fs= Ms.GetFittedStar();
  unsigned fsIndex = fs->IndexInMatrix();
  if (_fittingPos)
    {
      Indices.push_back(fsIndex);
      Indices.push_back(fsIndex+1);
    }
  // For securing the outlier removal, the next block is just useless 
  if (_fittingPM)
     {
       for (unsigned k=0; k<NPAR_PM; ++k) Indices.push_back(fsIndex+2+k);
     }
  /* Should not put the index of refaction stuff or we will not be
     able to remove more than 1 star at a time. */
}

//! contributions to derivatives of (presumambly) outlier terms. No discarding done.
void AstromFit::OutliersContributions(MeasuredStarList &MOutliers,
				      FittedStarList &FOutLiers,
				      TripletList &TList, 
				      Eigen::VectorXd &Grad)
{
  // contributions from measurement terms:
  for (auto i= MOutliers.begin(); i!= MOutliers.end(); ++i)
    {
      MeasuredStar &out = **i;
      MeasuredStarList tmp;
      tmp.push_back(&out);
      const CcdImage &ccd = *(out.ccdImage);
      LSDerivatives1(ccd, TList, Grad, &tmp);
    }
  LSDerivatives2(FOutLiers, TList, Grad);
}


//! Discards measurements and reference contributions contributing to the chi2 more than a cut, computed as <chi2>+NSigCut+rms(chi2) (statistics over contributions to the chi2). Returns the number of removed outliers. No refit done.
unsigned AstromFit::RemoveOutliers(const double &NSigCut,
				   const std::string &MeasOrRef)
{
  MeasuredStarList MSOutliers;
  FittedStarList FSOutliers;
  unsigned n = FindOutliers(NSigCut, MSOutliers, FSOutliers, MeasOrRef);
  RemoveMeasOutliers(MSOutliers);
  RemoveRefOutliers(FSOutliers);
  return n;
}
  


//! Find Measurements and references contributing more than a cut, computed as <chi2>+NSigCut+rms(chi2). The outliers are NOT removed, and no refit is done.
/*! After returning from here, there are still measurements that
  contribute above the cut, but their contribution should be
  evaluated after a refit before discarding them. */
unsigned AstromFit::FindOutliers(const double &NSigCut,
				 MeasuredStarList &MSOutliers,
				 FittedStarList &FSOutliers,
				 const std::string &MeasOrRef) const
{
  bool searchMeas = (MeasOrRef.find("Meas") != std::string::npos);
  bool searchRef = (MeasOrRef.find("Ref") != std::string::npos);

  CcdImageList &L=_assoc.ccdImageList;
  // collect chi2 contributions
  Chi2Vect chi2s;
  chi2s.reserve(_nMeasuredStars+_assoc.refStarList.size());
  // contributions from measurement terms:
  if (searchMeas)
    AccumulateStatImageList(_assoc.ccdImageList, chi2s);
  // and from reference terms
  if (searchRef)
    AccumulateStatRefStars(chi2s);

  // do some stat
  unsigned nval = chi2s.size();
  if (nval==0) return 0;
  sort(chi2s.begin(), chi2s.end());
  double median = (nval & 1)? chi2s[nval/2].chi2 :  
    0.5*(chi2s[nval/2-1].chi2 + chi2s[nval/2].chi2);
  // some more stats. should go into the class if recycled anywhere else
  double sum=0; double sum2 = 0;
  for (auto i=chi2s.begin(); i!=chi2s.end(); ++i) 
    {sum+= i->chi2;sum2+= sqr(i->chi2);}
  double average = sum/nval;
  double sigma = sqrt(sum2/nval - sqr(average));
  cout << "INFO : RemoveOutliers chi2 stat: mean/median/sigma " 
       << average << '/'<< median << '/' << sigma << endl;
  double cut = average+NSigCut*sigma;
  /* For each of the parameters, we will not remove more than 1
     measurement that contributes to constraining it. Keep track using
     of what we are touching using an integer vector. This is the
     trick that Marc Betoule came up to for outlier removals in "star
     flats" fits. */
  Eigen::VectorXi affectedParams(_nParTot);
  affectedParams.setZero();

  unsigned nOutliers = 0; // returned to the caller
  // start from the strongest outliers.
  for (auto i = chi2s.rbegin(); i != chi2s.rend(); ++i)
    {
      if (i->chi2 < cut) break; // because the array is sorted. 
      vector<unsigned> indices;
      indices.reserve(100); // just there to limit reallocations.
      /* now, we want to get the indices of the parameters this chi2
	 term depends on. We have to figure out which kind of term it
	 is; we use for that the type of the star attached to
	 the Chi2Entry. */
      MeasuredStar *ms = dynamic_cast<MeasuredStar *>(i->ps); 
      FittedStar *fs = NULL;
      if (!ms) // it is reference term.
	{
	  fs = dynamic_cast<FittedStar *>(i->ps); 
	  indices.push_back(fs->IndexInMatrix());
	  indices.push_back(fs->IndexInMatrix()+1); // probably useless
	  /* One might think it would be useful to account for PM
	     parameters here, but it is just useless */
	}
      else // it is a measurement term.
	{
	  GetMeasuredStarIndices(*ms, indices);
	}

      /* Find out if we already discarded a stronger outlier
	 constraining some parameter this one constrains as well. If
	 yes, we keep this one, because this stronger outlier could be 
	 causing the large chi2 we have in hand.  */
      bool drop_it = true;
      for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	if (affectedParams(*i) !=0) drop_it = false;
      
      if (drop_it) // store the outlier in one of the lists:
	{
	  if (ms) // measurement term
	    MSOutliers.push_back(ms);
	  else // ref term
	    FSOutliers.push_back(fs);
	  // mark the parameters as directly changed when we discard this chi2 term.
	  for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	    affectedParams(*i)++;
	  nOutliers++;
	}
    } // end loop on measurements/references
  cout << "INFO : FindOutliers : found " 
       << MSOutliers.size() << " meas outliers and " 
       << FSOutliers.size ()<< " ref outliers " << endl;
  
  return nOutliers;
}


void AstromFit::RemoveMeasOutliers(MeasuredStarList &Outliers)
{
  for (auto i = Outliers.begin(); i!= Outliers.end(); ++i)
    {
      MeasuredStar &ms = **i;
      FittedStar *fs = const_cast<FittedStar *>(ms.GetFittedStar());
      ms.SetValid(false); 
      fs->MeasurementCount()--; // could be put in SetValid
    }
}
  

void AstromFit::RemoveRefOutliers(FittedStarList &Outliers)
{
  for (auto i=Outliers.begin(); i!= Outliers.end() ;++i)
    {
      FittedStar &fs = **i;
      fs.SetRefStar(NULL);
    }
}



/*! WhatToFit is searched for strings : "Distortions", "Positions",
"Refrac", "PM" which define which parameter set is going to be
variable when computing derivatives (LSDerivatives) and minimizing
(Minimize()).  WhatToFit="Positions Distortions" will minimize w.r.t
mappings and objects positions, and not w.r.t proper motions and
refraction modeling.  However if proper motions and/or refraction
parameters have already been set, then they are accounted for when
computing residuals.  The string is forwarded to the DistortionModel,
and it can then be used to turn subsets of distortion parameter on or
off, if the DistortionModel implements such a thing.
*/
void AstromFit::AssignIndices(const std::string &WhatToFit)
{
  _WhatToFit = WhatToFit;
  cout << "INFO: we are going to fit : " << WhatToFit << endl;
  _fittingDistortions = (_WhatToFit.find("Distortions") != string::npos);
  _fittingPos = (_WhatToFit.find("Positions") != string::npos);
  _fittingRefrac = (_WhatToFit.find("Refrac") != string::npos);
  if (_sigCol == 0 && _fittingRefrac) 
    {
      cout << "WARNING: We cannot fit refraction coefficients without a color lever arm. Ignoring refraction" << endl;
      _fittingRefrac = false;
    }	
  _fittingPM = (_WhatToFit.find("PM") != string::npos);
// When entering here, we assume that WhatToFit has already been interpreted.


  _nParDistortions = 0;
  if (_fittingDistortions) 
    _nParDistortions = _distortionModel->AssignIndices(0,_WhatToFit);
  unsigned ipar = _nParDistortions;

  if (_fittingPos)
    {
      FittedStarList &fsl = _assoc.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when updating (OffsetParams())
	  // - in GetMeasuredStarIndices
	  fs.SetIndexInMatrix(ipar);
	  ipar+=2;
	  if ((_fittingPM) & fs.mightMove) ipar+= NPAR_PM;
	}
    }
  unsigned _nParPositions = ipar-_nParDistortions;
  if (_fittingRefrac)
    { 
      _refracPosInMatrix = ipar;
      ipar += _nRefrac;
    }
  _nParTot = ipar;

#if (0)  
  //DEBUG
  cout << " INFO: np(d,p, total) = " 
       <<  _nParDistortions << ' '
       << _nParPositions << ' '  
       << _nParTot << ' '
       << WhatToFit << endl;
  const FittedStar &ffs = (**(_assoc.fittedStarList.begin()));
  cout << " INFO : first Star Index : " <<  ffs.IndexInMatrix() << ' ' << Point(ffs) << endl;
#endif
}


      
void AstromFit::OffsetParams(const Eigen::VectorXd& Delta)
{
  if (Delta.size() != _nParTot) 
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "AstromFit::OffsetParams : the provided vector length is not compatible with the current WhatToFit setting");
  if (_fittingDistortions) 
    _distortionModel->OffsetParams(Delta);

  if (_fittingPos)
    {
      FittedStarList &fsl = _assoc.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when assigning indices (AssignIndices())
	  unsigned index = fs.IndexInMatrix();
	  fs.x += Delta(index);
	  fs.y += Delta(index+1);
	  if ((_fittingPM) & fs.mightMove)
	    {
	      fs.pmx += Delta(index+2);
	      fs.pmy += Delta(index+3);
	    }
	}
    }
  if (_fittingRefrac)
    {
      for (unsigned k=0; k<_nRefrac; ++k)
	_refracCoefficient[k] += Delta(_refracPosInMatrix+k);
    }
}

// should not be too large !
#ifdef STORAGE
static void write_sparse_matrix_in_fits(const SpMat &mat, const string &FitsName)
{
  if (mat.rows()*mat.cols() > 2e8)
    {
      cout << "WARNING :  write_sparse_matrix_in_fits : yout matrix is too large. " << FitsName << " not generated"<< endl;
      return;
    }
  Mat m(mat.rows(),mat.cols());
  for (int k=0; k<mat.outerSize(); ++k)
    for (SpMat::InnerIterator it(mat,k); it; ++it)
      {
	m (it.row(), it.col()) = it.value();
      }
  m.writeFits(FitsName);
}

static void write_vect_in_fits(const Eigen::VectorXd &V, const string &FitsName)
{
  Vect v(V.size());
  for (int k=0; k <V.size(); ++k) v(k) = V(k);
  Mat(v).writeFits(FitsName);
}

#endif


/*! This is a complete Newton Raphson step. Compute first and 
  second derivatives, solve for the step and apply it, without 
  a line search. */
unsigned AstromFit::Minimize(const std::string &WhatToFit, const double NSigRejCut)
{
  AssignIndices(WhatToFit);
  
  // return code can take 3 values : 
  // 0 : fit has converged - no more outliers
  // 1 : still some ouliers but chi2 increases
  // 2 : factorization failed
  unsigned returnCode = 0;

  // TODO : write a guesser for the number of triplets
  unsigned nTrip = (_LastNTrip) ? _LastNTrip: 1e6;
  TripletList tList(nTrip);
  Eigen::VectorXd grad(_nParTot);  grad.setZero();

  //Fill the triplets
  clock_t tstart = clock();
  LSDerivatives(tList, grad);
  clock_t tend = clock();
  _LastNTrip = tList.size(); 

  cout << " INFO: End of triplet filling, ntrip = " << tList.size() 
       << " CPU = " << float(tend-tstart)/float(CLOCKS_PER_SEC) 
       << endl;

  SpMat hessian;
  {
#if (TRIPLET_INTERNAL_COORD == COL)
    SpMat jacobian(_nParTot,tList.NextFreeIndex());
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    tList.clear(); //tList.shrink_to_fit();
    clock_t tstart = clock();
    hessian = jacobian*jacobian.transpose();
    clock_t tend = clock();
    std::cout << "INFO: CPU for J*Jt " 
	      << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;

#else
    SpMat jacobian(tList.NextRank(), _nParTot);
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    tList.clear(); //tList.shrink_to_fit(); 
    cout << " starting H=JtJ " << endl;
    hessian = jacobian.transpose()*jacobian;
#endif
  }// release the Jacobian


  cout << "INFO: hessian : dim=" << hessian.rows() 
       << " nnz=" << hessian.nonZeros() 
       << " filling-frac = " << hessian.nonZeros()/sqr(hessian.rows()) << endl;
  cout << "INFO: starting factorization" << endl;

  tstart = clock();
  CholmodSimplicialLDLT2<SpMat> chol(hessian);
  if (chol.info() != Eigen::Success)
    {
      cout << "ERROR: AstromFit::Minimize : factorization failed " << endl;
      return 2;
    }

  tend = clock();
  std::cout << "INFO: CPU for factorize-solve " 
  	    << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;
  tstart = tend;

  unsigned tot_outliers = 0;
  double old_chi2 = ComputeChi2().chi2;

  while (true)
    {
      Eigen::VectorXd delta = chol.solve(grad);
      //  cout << " offsetting parameters" << endl;
      OffsetParams(delta);
      Chi2 current_chi2(ComputeChi2());
      cout << current_chi2 << endl;
      if (current_chi2.chi2 > old_chi2)
	{
	  cout << "WARNING: chi2 went up, exiting outlier rejection loop" << endl;
	  returnCode = 1;
	  break;
	}
      old_chi2 = current_chi2.chi2;

      if (NSigRejCut == 0) break;
      MeasuredStarList moutliers;
      FittedStarList foutliers;
      int n_outliers = FindOutliers(NSigRejCut, moutliers, foutliers);
      tot_outliers += n_outliers;
      if (n_outliers == 0) break;
      TripletList tList(1000); // initial allocation size.
      grad.setZero(); // recycle the gradient
      // compute the contributions of outliers to derivatives
      OutliersContributions(moutliers, foutliers, tList, grad);
      // actually discard them
      RemoveMeasOutliers(moutliers);
      RemoveRefOutliers(foutliers);      
      // convert triplet list to eigen internal format
      SpMat h(_nParTot,tList.NextFreeIndex());
      h.setFromTriplets(tList.begin(), tList.end());
      int update_status = chol.update(h, false /* means downdate */);
      cout << "INFO: cholmod  update_status " << update_status << endl;
      /* The contribution of outliers to the gradient is the opposite
	 of the contribution of all other terms, because they add up
	 to 0 */
      grad *= -1;
      tend = clock();
      std::cout << "INFO: CPU for chi2-update_factor "  
		<< float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;
      tstart = tend;
    }

  cout << "INFO: total number of outliers " << tot_outliers << endl;

  return returnCode;
}

/* DEBUGGING routine */
void AstromFit::CheckStuff()
{
#if (0)
  const char *what2fit[] = {"Positions", "Distortions", "Refrac",
		      "Positions Distortions", "Positions Refrac", 
		      "Distortions Refrac",
		      "Positions Distortions Refrac"};
#endif
  const char *what2fit[] = {"Positions", "Distortions",
		      "Positions Distortions"};
  // DEBUG
  for (int k=0; k < sizeof(what2fit)/sizeof(what2fit[0]); ++k)
    {
      AssignIndices(what2fit[k]);
#if (0)
      fsIndexDebug = _nParDistortions;
      heavyDebug = true;
#endif
      TripletList tList(10000);
      Eigen::VectorXd rhs(_nParTot);  rhs.setZero();		    
      LSDerivatives(tList, rhs);
      SpMat jacobian(_nParTot,tList.NextFreeIndex());
      jacobian.setFromTriplets(tList.begin(), tList.end());
      SpMat hessian = jacobian*jacobian.transpose();
#ifdef STORAGE      
      char name[24];
      sprintf(name,"h%d.fits", k);
      write_sparse_matrix_in_fits(hessian, name);
      sprintf(name,"g%d.fits", k);
      write_vect_in_fits(rhs, name);
#endif
      cout << "npar : " << _nParTot << ' ' << _nParDistortions << ' ' << endl;

    }
}

void AstromFit::MakeResTuple(const std::string &TupleName) const
{
  /* cook-up 2 different file names by inserting something just before
     the dot (if any), and within the actual file name. */
  size_t dot = TupleName.rfind('.');
  size_t slash = TupleName.rfind('/');
  if (dot == string::npos || (slash != string::npos && dot < slash))
    dot = TupleName.size();  
  std::string meas_tuple(TupleName);
  meas_tuple.insert(dot,"-meas");
  MakeMeasResTuple(meas_tuple);
  std::string ref_tuple(TupleName);
  ref_tuple.insert(dot,"-ref");
  MakeRefResTuple(ref_tuple);
}

void AstromFit::MakeMeasResTuple(const std::string &TupleName) const
{
  std::ofstream tuple(TupleName.c_str());
  tuple << "#xccd: coordinate in CCD" << endl
	<< "#yccd: " << endl
	<< "#rx:   residual in degrees in TP" << endl
	<< "#ry:" << endl
	<< "#xtp: transformed coordinate in TP " << endl
	<< "#ytp:" << endl 
	<< "#mag: rough mag" << endl
	<< "#jd: Julian date of the measurement" << endl
    	<< "#rvx: transformed measurement uncertainty " << endl
    	<< "#rvy:" << endl
    	<< "#rvxy:" << endl
	<< "#color : " << endl
	<< "#fsindex: some unique index of the object" << endl
	<< "#ra: pos of fitted star" << endl
	<< "#dec: pos of fitted star" << endl
	<< "#chi2: contribution to Chi2 (2D dofs)" << endl
	<< "#nm: number of measurements of this FittedStar" << endl
    	<< "#chip: chip number" << endl
    	<< "#shoot: shoot id" << endl
	<< "#end" << endl;
  const CcdImageList &L=_assoc.TheCcdImageList();
  for (auto i=L.cbegin(); i!=L.end() ; ++i)
    {
      const CcdImage &im = **i;
      const MeasuredStarList &cat = im.CatalogForFit();
      const Mapping *mapping = _distortionModel->GetMapping(im);
      const Point &refractionVector = im.ParallacticVector();
      double jd = im.JD() - _JDRef;
      double zp = im.ZP();
      unsigned iband= im.BandRank();
      for (auto is=cat.cbegin(); is!=cat.end(); ++is)
	{
	  const MeasuredStar &ms = **is;
	  if (!ms.IsValid()) continue;
	  FatPoint tpPos;
	  FatPoint inPos = ms;
	  TweakAstromMeasurementErrors(inPos, ms, _posError);
	  mapping->TransformPosAndErrors(inPos, tpPos);
	  const Gtransfo* sky2TP = _distortionModel->Sky2TP(im);
	  const FittedStar *fs = ms.GetFittedStar();
	  
	  Point fittedStarInTP = TransformFittedStar(*fs, sky2TP,
						     refractionVector, 
						     _refracCoefficient[iband],
						     jd);
	  Point res=tpPos-fittedStarInTP;
	  double det = tpPos.vx*tpPos.vy-sqr(tpPos.vxy);
	  double wxx = tpPos.vy/det;
	  double wyy = tpPos.vx/det;
	  double wxy = -tpPos.vxy/det;
	  //	  double chi2 = rx*(wxx*rx+wxy*ry)+ry*(wxy*rx+wyy*ry);
	  double chi2 = wxx*res.x*res.x + wyy*res.y*res.y + 2*wxy*res.x*res.y;
	  tuple << std::setprecision(9);
	  tuple << ms.x << ' ' << ms.y << ' ' 
		<< res.x << ' ' << res.y << ' '
		<< tpPos.x << ' ' << tpPos.y << ' '
		<< fs->Mag() << ' ' << jd << ' ' 
		<< tpPos.vx << ' ' << tpPos.vy << ' ' << tpPos.vxy << ' ' 
		<< fs->color << ' ' 
		<< fs->IndexInMatrix() << ' '
	        << fs->x << ' ' << fs->y << ' '
		<< chi2 << ' ' 
		<< fs->MeasurementCount() << ' ' 
		<< im.Chip() << ' ' << im.Shoot() << endl;
	}// loop on measurements in image
    }// loop on images

}

void AstromFit::MakeRefResTuple(const std::string &TupleName) const
{
  std::ofstream tuple(TupleName.c_str());
  tuple << "#ra: coordinates of FittedStar" << endl
	<< "#dec: " << endl
	<< "#rx:   residual in degrees in TP" << endl
	<< "#ry:" << endl
	<< "#mag: mag" << endl
    	<< "#rvx: transformed measurement uncertainty " << endl
    	<< "#rvy:" << endl
    	<< "#rvxy:" << endl
	<< "#color : " << endl
	<< "#fsindex: some unique index of the object" << endl
	<< "#chi2: contribution to Chi2 (2D dofs)" << endl
	<< "#nm: number of measurements of this FittedStar" << endl
	<< "#end" << endl;
  // The following loop is heavily inspired from AstromFit::ComputeChi2()
  const FittedStarList &fsl = _assoc.fittedStarList;
  TanRaDec2Pix proj(GtransfoLin(), Point(0.,0.));
  for (auto i = fsl.cbegin(); i!= fsl.end(); ++i)
    {
      const FittedStar &fs = **i;
      const RefStar *rs = fs.GetRefStar();
      if (rs == NULL) continue;
      proj.SetTangentPoint(fs);
      // fs projects to (0,0), no need to compute its transform.
      FatPoint rsProj;
      proj.TransformPosAndErrors(*rs, rsProj);
      double rx = rsProj.x; // -fsProj.x (which is 0)
      double ry = rsProj.y;
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      double wxx = rsProj.vy/det;
      double wyy = rsProj.vx/det;
      double wxy = -rsProj.vxy/det;
      double chi2 = wxx*sqr(rx) + 2*wxy*rx*ry+ wyy*sqr(ry);
      tuple << std::setprecision(9);
      tuple << fs.x << ' ' << fs.y << ' ' 
	    << rx << ' ' << ry << ' '
	    << fs.Mag() << ' ' 
	    << rsProj.vx << ' ' << rsProj.vy << ' ' << rsProj.vxy << ' ' 
	    << fs.color << ' ' 
	    << fs.IndexInMatrix() << ' '
	    << chi2 << ' ' 
	    << fs.MeasurementCount() << endl;
    }// loop on FittedStars
}


}} // end of namespaces
