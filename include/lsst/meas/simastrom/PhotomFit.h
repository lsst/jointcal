#ifndef PHOTOMFIT__H
#define PHOTOMFIT__H

#include <string>
#include <iostream>
#include <sstream>
#include <map>

#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Eigenstuff.h"
#include "lsst/meas/simastrom/Tripletlist.h"
#include "lsst/meas/simastrom/PhotomModel.h"
#include "lsst/meas/simastrom/Chi2.h"

namespace lsst {
namespace meas {
namespace simastrom {

class Associations;

/*! Some comments.
  

*/


//! Class that handles the photometric least squares problem.

class PhotomFit {
  private : 

  Associations &_assoc; 
  std::string _WhatToFit;
  bool _fittingModel, _fittingFluxes;
  unsigned _nParModel, _nParFluxes, _nParTot;
  PhotomModel * _photomModel;
  double _fluxError;
  int _LastNTrip; // last triplet count, used to speed up allocation


  
 public :

  //! this is the only constructor
  PhotomFit (Associations &A, PhotomModel *M, double PosError);
  
  //! Does a 1 step minimization, assuming a linear model.
  /*! It calls AssignIndices, LSDerivatives, solves the linear system
    and calls OffsetParams. No line search. Relies on sparse linear
    algebra. */
  bool Minimize(const std::string &WhatToFit);

  //! Derivatives of the Chi2
  void LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Compute the derivatives for this CcdImage
  void LSDerivatives(const CcdImage &Ccd, 
		     TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Set parameter groups fixed or variable and assign indices to each parameter in the big matrix (which will be used by OffsetParams(...).
  void AssignIndices(const std::string &WhatToFit);

  //! Offsest the parameters by the requested quantities. The used parameter layout is the one from the last call to AssignIndices or Minimize().
  /*! There is no easy way to check that the current setting of
      WhatToFit and the provided Delta vector are compatible. We can
      only test the size. */
  void OffsetParams(const Eigen::VectorXd &Delta);

  //! Returns a chi2 for the current state
  Chi2 ComputeChi2() const;

  //! Produces an ntuple 
  void MakeResTuple(const std::string &TupleName) const;

 private:

  template <class ListType, class Accum> 
    void AccumulateStat(ListType &L, Accum &A) const;


#ifdef STORAGE
  //! Compute the derivatives of the reference terms
  void LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Evaluates the chI^2 derivatives (Jacobian and gradient) for the current WhatToFit setting.


  //! returns how many outliers were removed. No refit done.
  unsigned RemoveOutliers(const double &NSigCut);

  //! Produces a tuple containing residuals of measurement terms.
  void MakeMeasResTuple(const std::string &TupleName) const;

  //! Produces a tuple containing residuals of reference terms.
  void MakeRefResTuple(const std::string &TupleName) const;


  //! access to the fitted refraction coefficients. Unit depends on scale in the tangentPlane. Degrees for an actual tangent plane.
  std::vector<double> RefractionCoefficients() const 
    { return _refracCoefficient;}

  void CheckStuff();

 private : 

  Point TransformFittedStar(const FittedStar &F,
			    const Gtransfo * Sky2TP,
			    const Point &RefractionVector,
			    const double &RefractionCoeff,
			    const double &Jd) const;

  //! only for outlier removal
  void GetMeasuredStarIndices(const MeasuredStar &Ms, 
			      std::vector<unsigned> &Indices) const;
#endif
};


}}}
#endif /* PHOTOMFIT__H */
