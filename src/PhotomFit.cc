#include <iostream>
#include <iomanip>
#include <algorithm>
#include "lsst/jointcal/PhotomFit.h"
#include "lsst/jointcal/Associations.h"

#include "lsst/jointcal/Gtransfo.h"
#include "Eigen/Sparse"
//#include "Eigen/CholmodSupport" // to switch to cholmod
#include "lsst/pex/exceptions.h"
#include <fstream>
#include "lsst/jointcal/Tripletlist.h"

typedef Eigen::SparseMatrix<double> SpMat;

using namespace std;

static double sqr(double x) {return x*x;}

namespace lsst {
namespace jointcal {


PhotomFit::PhotomFit(Associations &associations, PhotomModel *photomModel, double fluxError) :
  _associations(associations),  _photomModel(photomModel), _fluxError(fluxError),
  _lastNTrip(0)
{
  // The various _npar... are initialized in assignIndices.
  // Although there is no reason to adress them before one might be tempted by
  // evaluating a Chi2 rightaway, .. which uses these counts, so:
  assignIndices("");
}

void PhotomFit::LSDerivatives(TripletList &tList, Eigen::VectorXd &rhs) const
{
  auto L = _associations.getCcdImageList();
  for (auto im=L.cbegin(); im!=L.end() ; ++im)
    {
      LSDerivatives(**im, tList, rhs);
    }
}

void PhotomFit::LSDerivatives(const CcdImage &ccdImage,
                              TripletList &tList, Eigen::VectorXd &rhs,
                              const MeasuredStarList *measuredStarList) const
{
  /***************************************************************************/
  /**  Changes in this routine should be reflected into accumulateStat       */
  /***************************************************************************/
  /* this routine works in two different ways: either providing the
     Ccd, of providing the MeasuredStarList. In the latter case, the
     Ccd should match the one(s) in the list. */
  if (measuredStarList)
    assert ( (*(measuredStarList->begin()))->ccdImage == &ccdImage);

  unsigned npar_max = 100; // anything large
  vector<unsigned> indices(npar_max,-1);

  Eigen::VectorXd h(npar_max);
  Eigen::VectorXd grad(npar_max);
  // current position in the Jacobian
  unsigned kTriplets = tList.NextFreeIndex();
  const MeasuredStarList &catalog = (measuredStarList) ? *measuredStarList : ccdImage.getCatalogForFit();

  for (auto i = catalog.begin(); i!= catalog.end(); ++i)
    {
      const MeasuredStar& measuredStar = **i;
      if (!measuredStar.IsValid()) continue;
      // tweak the measurement errors
      double sigma=measuredStar.eflux;
#ifdef FUTURE
      TweakPhotomMeasurementErrors(inPos, measuredStar, _posError);
#endif
      h.setZero(); // we cannot be sure that all entries will be overwritten.

      double pf = _photomModel->photomFactor(ccdImage, measuredStar);
      const FittedStar *fs = measuredStar.GetFittedStar();

      double res = measuredStar.flux - pf * fs->flux;

      if (_fittingModel)
	{
	  _photomModel->getIndicesAndDerivatives(measuredStar,
                                             ccdImage,
                                             indices,
                                             h);
	  for (unsigned k=0; k<indices.size(); k++)
	    {
	      unsigned l = indices[k];
	      tList.AddTriplet(l, kTriplets, h[k]*fs->flux/sigma);
	      rhs[l] += h[k]*res/sqr(sigma);
	    }
	}
      if (_fittingFluxes)
	{
	  unsigned index = fs->IndexInMatrix();
	  tList.AddTriplet(index,kTriplets, pf/sigma);
	  rhs[index] += res*pf/sqr(sigma);
	}
      kTriplets += 1; // each measurement contributes 1 column in the Jacobian
    } // end loop on measurements
  tList.SetNextFreeIndex(kTriplets);
}

// This is almost a selection of lines of LSDerivatives(CcdImage ...)
/* This routine is template because it is used
both with its first argument as "const CcdImageList &" and "CcdImageList &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ListType, class Accum>
void PhotomFit::accumulateStat(ListType &listType, Accum &accum) const
{
  for (auto im=listType.begin(); im!=listType.end() ; ++im)
    {
  /**********************************************************************/
  /**  Changes in this routine should be reflected into LSDerivatives  */
  /**********************************************************************/
      auto &ccdIMage = **im;
      auto &catalog = ccdIMage.getCatalogForFit();

      for (auto i = catalog.begin(); i!= catalog.end(); ++i)
	{
	  auto &measuredStar = **i;
	  if (!measuredStar.IsValid()) continue;
	  // tweak the measurement errors
	  double sigma=measuredStar.eflux;
#ifdef FUTURE
	  TweakPhotomMeasurementErrors(inPos, measuredStar, _posError);
#endif

	  double pf = _photomModel->photomFactor(ccdIMage, measuredStar);
	  const FittedStar *fs = measuredStar.GetFittedStar();
	  double res = measuredStar.flux - pf * fs->flux;
	  double chi2Val = sqr(res/sigma);
	  accum.AddEntry(chi2Val, 1, &measuredStar);
	} // end loop on measurements
    }
}

//! for the list of images in the provided  association and the reference stars, if any
Chi2 PhotomFit::computeChi2() const
{
  Chi2 chi2;
  accumulateStat(_associations.getCcdImageList(), chi2);
  // so far, chi2.ndof contains the number of squares.
  // So, subtract here the number of parameters.
  chi2.ndof -= _nParTot;
  return chi2;
}

void PhotomFit::outliersContributions(MeasuredStarList &outliers,
                                      TripletList &tList,
                                      Eigen::VectorXd &grad)
{
  for (auto i= outliers.begin(); i!= outliers.end(); ++i)
    {
      MeasuredStar &out = **i;
      MeasuredStarList tmp;
      tmp.push_back(&out);
      const CcdImage &ccdImage = *(out.ccdImage);
      LSDerivatives(ccdImage, tList, grad, &tmp);
      out.SetValid(false);
      FittedStar *fs = const_cast<FittedStar *>(out.GetFittedStar());
      fs->MeasurementCount()--;
    }
}


//! a class to accumulate chi2 contributions together with pointers to the contributors.
/*! This structure allows to compute the chi2 statistics (average and
  variance) and directly point back to the bad guys without
  relooping. The Chi2Entry routine makes it compatible with
  accumulateStatImage and accumulateStatImageList. */
struct Chi2Entry
{
  double chi2;
  MeasuredStar *measuredStar;

  Chi2Entry(double c, MeasuredStar *s): chi2(c), measuredStar(s) {}
  // for sort
  bool operator < (const Chi2Entry &R) const {return (chi2<R.chi2);}
};

struct Chi2Vect : public vector<Chi2Entry>
{
  void AddEntry(double Chi2Val, unsigned ndof, MeasuredStar *measuredStar)
  { push_back(Chi2Entry(Chi2Val, measuredStar));}

};

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void PhotomFit::getMeasuredStarIndices(const MeasuredStar &measuredStar,
                                       std::vector<unsigned> &indices) const
{
  indices.clear();
  if (_fittingModel)
    {
      Eigen::VectorXd h(100);
      _photomModel->getIndicesAndDerivatives(measuredStar,
                                             *measuredStar.ccdImage,
                                             indices,
                                             h);
    }
  if (_fittingFluxes)
    {
      const FittedStar *fs= measuredStar.GetFittedStar();
      unsigned fsIndex = fs->IndexInMatrix();
      indices.push_back(fsIndex);
    }
}

void PhotomFit::findOutliers(double nSigCut, MeasuredStarList &outliers) const
{
  /* Aims at providing an outlier list for small-rank update
     of the factorization. */

  // collect chi2 contributions of the measurements.
  Chi2Vect chi2s;
  //  chi2s.reserve(_nMeasuredStars);
  accumulateStat(_associations.ccdImageList, chi2s);
  // do some stat
  unsigned nval = chi2s.size();
  if (nval==0) return;
  sort(chi2s.begin(), chi2s.end()); // increasing order. We rely on this later.
  double median = (nval & 1)? chi2s[nval/2].chi2 :
    0.5*(chi2s[nval/2-1].chi2 + chi2s[nval/2].chi2);
  // some more stats. should go into the class if recycled anywhere else
  double sum=0; double sum2 = 0;
  for (auto i=chi2s.begin(); i!=chi2s.end(); ++i)
    {sum+= i->chi2;sum2+= sqr(i->chi2);}
  double average = sum/nval;
  double sigma = sqrt(sum2/nval - sqr(average));
  cout << "INFO : findOutliers chi2 stat: mean/median/sigma "
       << average << '/'<< median << '/' << sigma << endl;
  double cut = average+nSigCut*sigma;
  /* For each of the parameters, we will not remove more than 1
     measurement that contributes to constraining it. Keep track
     of the affected parameters using an integer vector. This is the
     trick that Marc Betoule came up to for outlier removals in "star
     flats" fits. */
  Eigen::VectorXi affectedParams(_nParTot);
  affectedParams.setZero();

  // start from the strongest outliers, i.e. at the end of the array.
  for (auto i = chi2s.rbegin(); i != chi2s.rend(); ++i)
    {
      if (i->chi2 < cut) break; // because the array is sorted.
      vector<unsigned> indices;
      getMeasuredStarIndices(*(i->measuredStar), indices);
      bool drop_it = true;
      /* find out if a stronger outlier constraining one of the parameters
	 this one contrains was already discarded. If yes, we keep this one */
      for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	if (affectedParams(*i) !=0) drop_it = false;

      if (drop_it)
	{
	  for (auto i=indices.cbegin(); i!= indices.end(); ++i)
	    affectedParams(*i)++;
	  outliers.push_back(i->measuredStar);
	}
    } // end loop on measurements
  cout << "INFO : findMeasOutliers : found "
       << outliers.size() << " outliers" << endl;
}

void PhotomFit::assignIndices(const std::string &whatToFit)
{
  _whatToFit = whatToFit;
  cout << "INFO: we are going to fit : " << whatToFit << endl;
  _fittingModel = (_whatToFit.find("Model") != string::npos);
  _fittingFluxes = (_whatToFit.find("Fluxes") != string::npos);
// When entering here, we assume that whatToFit has already been interpreted.

  _nParModel =   (_fittingModel) ? _photomModel->assignIndices(whatToFit,0) : 0;
  unsigned ipar = _nParModel;

  if (_fittingFluxes)
    {
      FittedStarList &fsl = _associations.fittedStarList;
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

void PhotomFit::offsetParams(const Eigen::VectorXd& delta)
{
  if (delta.size() != _nParTot)
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "PhotomFit::OffsetParams : the provided vector length is not compatible with the current whatToFit setting");
  if (_fittingModel)
    _photomModel->offsetParams(delta);

  if (_fittingFluxes)
    {
      FittedStarList &fsl = _associations.fittedStarList;
      for (FittedStarIterator i= fsl.begin(); i != fsl.end(); ++i)
	{
	  FittedStar &fs = **i;
	  // the parameter layout here is used also
	  // - when filling the derivatives
	  // - when assigning indices (assignIndices())
	  unsigned index = fs.IndexInMatrix();
	  fs.flux += delta(index);
	}
    }
}

bool PhotomFit::minimize(const std::string &whatToFit)
{
  assignIndices(whatToFit);

  // TODO : write a guesser for the number of triplets
  unsigned nTrip = (_lastNTrip) ? _lastNTrip: 1e6;
  TripletList tList(nTrip);
  Eigen::VectorXd grad(_nParTot);  grad.setZero();

  //Fill the triplets
  LSDerivatives(tList, grad);
  _lastNTrip = tList.size();

  SpMat hessian;
  {
    SpMat jacobian(_nParTot,tList.NextFreeIndex());
    jacobian.setFromTriplets(tList.begin(), tList.end());
    // release memory shrink_to_fit is C++11
    tList.clear(); //tList.shrink_to_fit();
    hessian = jacobian*jacobian.transpose();
  }// release the Jacobian


  cout << "INFO: hessian : dim=" << hessian.rows()
       << " nnz=" << hessian.nonZeros()
       << " filling-frac = " << hessian.nonZeros()/sqr(hessian.rows()) << endl;
  cout << "INFO: starting factorization" << endl;

  Eigen::SimplicialLDLT<SpMat> chol(hessian);
  if (chol.info() != Eigen::Success)
    {
      cout << "ERROR: PhotomFit::minimize : factorization failed " << endl;
      return false;
    }

  Eigen::VectorXd delta = chol.solve(grad);

  offsetParams(delta);
  return true;
}


void PhotomFit::makeResTuple(const std::string &tupleName) const
{
  std::ofstream tuple(tupleName.c_str());
  /* If we think the some coordinate on the focal plane is relevant in
     the ntuple, because thmodel relies on it, then we have to add
     some function to the model that returns this relevant
     coordinate. */
  tuple << "#xccd: coordinate in CCD" << endl
	<< "#yccd: " << endl
	<< "#mag: rough mag" << endl
	<< "#flux : measured flux" << endl
	<< "#eflux : measured flux erro" << endl
	<< "#fflux : fitted flux" << endl
	<< "#phot_factor:" <<  endl
	<< "#jd: Julian date of the measurement" << endl
	<< "#color : " << endl
	<< "#fsindex: some unique index of the object" << endl
	<< "#ra: pos of fitted star" << endl
	<< "#dec: pos of fitted star" << endl
	<< "#chi2: contribution to Chi2 (1 dof)" << endl
	<< "#nm: number of measurements of this FittedStar" << endl
    	<< "#chip: chip number" << endl
    	<< "#visit: visit id" << endl
	<< "#end" << endl;
  const CcdImageList &L=_associations.getCcdImageList();
  for (auto i=L.cbegin(); i!=L.end() ; ++i)
    {
      const CcdImage &im = **i;
      const MeasuredStarList &cat = im.getCatalogForFit();
      for (auto is=cat.cbegin(); is!=cat.end(); ++is)
	{
	  const MeasuredStar &ms = **is;
	  if (!ms.IsValid()) continue;
	  double sigma = ms.eflux;
#ifdef FUTURE
	  tweakPhotomMeasurementErrors(inPos, ms, _posError);
#endif
	  double pf = _photomModel->photomFactor(im, ms);
	  double jd = im.getMjd();
	  const FittedStar *fs = ms.GetFittedStar();
	  double res = ms.flux - pf * fs->flux;
	  double chi2Val = sqr(res/sigma);
	  tuple << ms.x << ' ' << ms.y << ' '
		<< fs->Mag() << ' '
		<< ms.flux << ' ' << ms.eflux << ' '
		<< fs->flux << ' '
		<< pf << ' '
		<< jd << ' '
		<< fs->color << ' '
		<< fs->IndexInMatrix() << ' '
	        << fs->x << ' ' << fs->y << ' '
		<< chi2Val << ' '
		<< fs->MeasurementCount() << ' '
		<< im.getCcdId() << ' ' << im.getVisit() << endl;
	}// loop on measurements in image
    }// loop on images

}

}} // end of namespaces
