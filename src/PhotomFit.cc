#include <iostream>
#include <iomanip>
#include <algorithm>
#include "lsst/meas/simastrom/PhotomFit.h"
#include "lsst/meas/simastrom/Associations.h"

#include "lsst/meas/simastrom/Gtransfo.h"
#include "Eigen/Sparse"
//#include "Eigen/CholmodSupport" // to switch to cholmod
#include <time.h> // for clock
#include "lsst/pex/exceptions.h"
#include <fstream>
#include "lsst/meas/simastrom/Tripletlist.h"

typedef Eigen::SparseMatrix<double> SpMat;

using namespace std;

static double sqr(const double &x) {return x*x;}

namespace lsst {
namespace meas {
namespace simastrom {


PhotomFit::PhotomFit(Associations &A, PhotomModel *M, double FluxError) : 
  _assoc(A),  _photomModel(M), _fluxError(FluxError)
{
  _LastNTrip = 0;

  //  _posError = PosError;

  //  _nMeasuredStars = 0;
  // The various _npar... are initialized in AssignIndices.
  // Although there is no reason to adress them before one might be tempted by
  // evaluating a Chi2 rightaway, .. which uses these counts, so:
  //  AssignIndices("");

}




/*! this is the first implementation of an error "model". 
  We'll certainly have to upgrade it. MeasuredStar provided
in case we need the mag.  */

//  static double posErrorIncrement = 0.02; // pixels
//  static double posErrorIncrement = 0.00; // pixels

#ifdef STORAGE
static void TweakPhotomMeasurementErrors(const MeasuredStar &Ms, double error)
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
#endif


void PhotomFit::LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs) const
{
  auto L = _assoc.TheCcdImageList();
  for (auto im=L.cbegin(); im!=L.end() ; ++im)
    {
      LSDerivatives(**im, TList, Rhs);
    }
}



// we could consider computing the chi2 here.
// (although it is not extremely useful)
void PhotomFit::LSDerivatives(const CcdImage &Ccd, 
			      TripletList &TList, Eigen::VectorXd &Rhs) const
{
  /***************************************************************************/
  /**  Changes in this routine should be reflected into AccumulateStatImage  */
  /***************************************************************************/
  /* Setup */
  // count parameters
  //  unsigned npar_flux = (_fittingFlux) ? 2 : 0;
  //  unsigned npar_refrac = (_fittingRefrac) ? 1 : 0;
  //  unsigned npar_pm = (_fittingPM) ? NPAR_PM : 0;
  //  unsigned npar_tot =  npar_mapping + npar_pos + npar_refrac + npar_pm;
  // if (npar_tot == 0) this CcdImage does not contribute 
  // any constraint to the fit, so :
  //  if (npar_tot == 0) return;
  unsigned npar_max = 100;
  vector<unsigned> indices(npar_max,-1);

  // the shape of h (et al) is required this way in order to be able to 
  // separate derivatives along x and y as vectors.
  Eigen::VectorXd h(npar_max);// halpha(npar_tot,2), hw(npar_tot,2); 
  Eigen::VectorXd grad(npar_max);
  // current position in the Jacobian
  unsigned kTriplets = TList.NextFreeIndex();
  const MeasuredStarList &catalog = Ccd.CatalogForFit();

  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      const MeasuredStar& ms = **i;
      if (!ms.IsValid()) continue;
      // tweak the measurement errors
      double sigma=ms.eflux;
#ifdef FUTURE
      TweakPhotomMeasurementErrors(inPos, ms, _posError);
#endif
      h.setZero(); // we cannot be sure that all entries will be overwritten.

      double pf = _photomModel->PhotomFactor(ms,Ccd);
      const FittedStar *fs = ms.GetFittedStar();

      double res = ms.flux - pf * fs->flux;
            
      if (_fittingModel) 
	{
	  unsigned nIndices = _photomModel->GetIndicesAndDerivatives(ms,
								     Ccd, 
								     indices,
								     h);
	  for (unsigned k=0; k<nIndices; k++)
	    {
	      unsigned l = indices[k];
	      TList.AddTriplet(l, kTriplets, h[k]*fs->flux/sigma);
	      Rhs[l] += h[k]*res/sqr(sigma);
	    }
	}
      if (_fittingFluxes)
	{
	  unsigned index = fs->IndexInMatrix();
	  TList.AddTriplet(index,kTriplets, pf/sigma);
	  Rhs[index] += res*pf/sqr(sigma);
	}
      kTriplets += 1; // each measurement contributes 1 column in the Jacobian
    } // end loop on measurements
  TList.SetNextFreeIndex(kTriplets);
}    

#ifdef STORAGE

// This is almost a selection of lines of LSDerivatives1(CcdImage ...)
/* This routine (and the following one) is template because it is used
both with its first argument as "const CCdImage &" and "CcdImage &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ImType, class Accum> 
void PhotomFit::AccumulateStatImage(ImType &Ccd, Accum &Accu) const
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
      TweakPhotomMeasurementErrors(inPos, ms, _posError);

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
void PhotomFit::AccumulateStatImageList(ListType &L, Accum &Accu) const
{
  for (auto im=L.begin(); im!=L.end() ; ++im)
    {
      AccumulateStatImage(**im, Accu);
    }
}

//! for the list of images in the provided  association and the reference stars, if any
Chi2 PhotomFit::ComputeChi2() const
{
  Chi2 chi2;
  AccumulateStatImageList(_assoc.TheCcdImageList(), chi2);
  // Now add the ref stars
  const FittedStarList &fsl = _assoc.fittedStarList;
  /* If you wonder why we project here, read comments in 
     PhotomFit::LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) */
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
      // TO DO : account for proper motions. If you do it here, think about MakeRefResTuple.
      double rx = rsProj.x; // -fsProj.x (which is 0)
      double ry = rsProj.y;
      double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
      double wxx = rsProj.vy/det;
      double wyy = rsProj.vx/det;
      double wxy = -rsProj.vxy/det;
      wxx *= HACK_REF_ERRORS;
      wyy *= HACK_REF_ERRORS;
      wxy *= HACK_REF_ERRORS;
      chi2.AddEntry(wxx*sqr(rx) + 2*wxy*rx*ry+ wyy*sqr(ry), 2, NULL);
    }
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
  MeasuredStar *ms;

  Chi2Entry(const double &c, MeasuredStar *s): chi2(c), ms(s) {}
  // for sort
  bool operator < (const Chi2Entry &R) const {return (chi2<R.chi2);}
};

struct Chi2Vect : public vector<Chi2Entry>
{
  void AddEntry(const double &Chi2Val, unsigned ndof, MeasuredStar *ms)
  { push_back(Chi2Entry(Chi2Val,ms));}

};
#endif


#ifdef STORAGE
unsigned PhotomFit::RemoveOutliers(const double &NSigCut)
{
  /* Some reshuffling would be needed if we wish to use the small-rank
     update trick rather than solving again. Typically We would
     need to compute the Jacobian and RHS contributions of the
     discarded measurement and update the current factorization and
     solution. */
  CcdImageList &L=_assoc.ccdImageList;
  // collect chi2 contributions
  Chi2Vect chi2s;
  chi2s.reserve(_nMeasuredStars);
  AccumulateStatImageList(_assoc.ccdImageList, chi2s);
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

  unsigned removed = 0; // returned to the caller
  // start from the strongest outliers.
  for (auto i = chi2s.rbegin(); i != chi2s.rend(); ++i)
    {
      if (i->chi2 < cut) break; // because the array is sorted. 
      vector<unsigned> indices;
      GetMeasuredStarIndices(*(i->ms), indices);
      bool drop_it = true;
      /* find out if a stronger outlier constraining one of the parameters
	 this one contrains was already discarded. If yes, we keep this one */
      for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	if (affectedParams(*i) !=0) drop_it = false;
      
      if (drop_it)
	{
	  FittedStar *fs = i->ms->GetFittedStar();
	  i->ms->SetValid(false); removed++;
	  fs->MeasurementCount()--; // could be put in SetValid
	  /* By making sure that we do not remove all MeasuredStars
	     pointing to a FittedStar in a single go,
	     fs->MeasurementCount() should never go to 0. 
	     
	     It seems plausible that the adopted mechanism prevents as
	     well to end up with under-constrained transfos. */
	  for (auto i=indices.cbegin(); i!= indices.cend(); ++i)
	    affectedParams(*i)++;
	}
    } // end loop on measurements
  cout << "INFO : RemoveOutliers : found and removed " 
       << removed << " outliers" << endl;
  return removed;
}
#endif

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
void PhotomFit::AssignIndices(const std::string &WhatToFit)
{
  _WhatToFit = WhatToFit;
  cout << "INFO: we are going to fit : " << WhatToFit << endl;
  _fittingModel = (_WhatToFit.find("Model") != string::npos);
  _fittingFluxes = (_WhatToFit.find("Fluxes") != string::npos);
// When entering here, we assume that WhatToFit has already been interpreted.


  if (_fittingModel) 
    _nParModel = _photomModel->AssignIndices(0);
  unsigned ipar = _nParModel;

  if (_fittingFluxes)
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
	  ipar+=1;
	}
    }
  _nParTot = ipar;
}

void PhotomFit::OffsetParams(const Eigen::VectorXd& Delta)
{
  if (Delta.size() != _nParTot) 
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "PhotomFit::OffsetParams : the provided vector length is not compatible with the current WhatToFit setting");
  if (_fittingModel) 
    _photomModel->OffsetParams(Delta);

  if (_fittingFluxes)
    {
      FittedStarList &fsl = _assoc.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when assigning indices (AssignIndices())
	  unsigned index = fs.IndexInMatrix();
	  fs.flux += Delta(index);
	}
    }
}


/*! This is a complete Newton Raphson step. Compute first and 
  second derivatives, solve for the step and apply it, without 
  a line search. */
bool PhotomFit::Minimize(const std::string &WhatToFit)
{
  AssignIndices(WhatToFit);

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
    SpMat jacobian(_nParTot,tList.NextFreeIndex());
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    tList.clear(); //tList.shrink_to_fit();
    clock_t tstart = clock();
    hessian = jacobian*jacobian.transpose();
    clock_t tend = clock();
    std::cout << "INFO: CPU for J*Jt " 
	      << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;
  }// release the Jacobian


  cout << "INFO: hessian : dim=" << hessian.rows() 
       << " nnz=" << hessian.nonZeros() 
       << " filling-frac = " << hessian.nonZeros()/sqr(hessian.rows()) << endl;
  cout << "INFO: starting factorization" << endl;

  tstart = clock();
  Eigen::SimplicialLDLT<SpMat> chol(hessian);
  if (chol.info() != Eigen::Success)
    {
      cout << "ERROR: PhotomFit::Minimize : factorization failed " << endl;
      return false;
    }

  Eigen::VectorXd delta = chol.solve(grad);

  //  cout << " offsetting parameters" << endl;
  OffsetParams(delta);
  tend = clock();
  std::cout << "INFO: CPU for factor-solve-update " 
  	    << float(tend-tstart)/float(CLOCKS_PER_SEC) << std::endl;
  return true;
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




/* DEBUGGING routine */
void PhotomFit::CheckStuff()
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

void PhotomFit::MakeResTuple(const std::string &TupleName) const
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

void PhotomFit::MakeMeasResTuple(const std::string &TupleName) const
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
	  TweakPhotomMeasurementErrors(inPos, ms, _posError);
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

void PhotomFit::MakeRefResTuple(const std::string &TupleName) const
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
  // The following loop is heavily inspired from PhotomFit::ComputeChi2()
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
	    << fs.MeasurementCount() << ' '
	    << endl;
    }// loop on FittedStars
}
#endif

}}}
