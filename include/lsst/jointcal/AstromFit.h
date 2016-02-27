#ifndef ASTROMFIT__H
#define ASTROMFIT__H

#include <string>
#include <iostream>
#include <sstream>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/Tripletlist.h"
#include "lsst/jointcal/DistortionModel.h"
#include "lsst/jointcal/Chi2.h"

namespace lsst {
namespace jointcal {

class Associations;


/*! This is the class that actually computes the quantities required
to carry out a LS astrometric fit wrt distortion mappings and coordinates
of common objects.  Namely it computes the Jacobian and
gradient of the chi2 (w.r.t. parameters), and the Chi2 itself.
It interfaces with the actual modelling of distortions 
via a mimimum virtual interface DistortionModel, and the actual mappings
via an other virtual interface : Mapping. 

  In short AstromFit aims at computing derivatives of 
least quares. The terms of the chi2 are of two kinds:

kind 1 ->   (T(X_M) - p(F))^T W (T(X_M) - p(F))

with X_M is a measured (2d) position in CCD coordinates, F refers to
the position of the object in some space, defined in practise by p.
There is one such term per measurement. The default setup would be
that p is the projection from sky to some tangent plane and hence T
maps the CCD coordinates onto this TP. p is obtained via the DistorsionModel
and can be different for all CcdImage's. Depending on what is beeing fitted,
one could imagine cases where the projector p is the same for all CcdImages.

Kind 2  -> (p'(F)-p'(R))^T  W_R (p'(F)-p'(R)) 
R refers to some  externally-provided reference object position, 
and p' to some projector from sky to some plane. The reference 
objects define the overall coordinate frame, which is required 
when all T and all F are fitted  simultaneously. There is one 
such term per external reference object. There can be more 
F (fitted) objects than R (reference) objects.

In the same framework, one can fit relative transforms between images by 
setting p = Identity for all input CcdImages and not fitting T for
one of the CcdImage's. One does not need reference object and 
would then naturally not have any Kind 2 terms.

*/



//! Class that handles the astrometric least squares problem.

class AstromFit {
  private : 

  Associations &_assoc; 
  std::string _WhatToFit;
  bool _fittingDistortions, _fittingPos, _fittingRefrac, _fittingPM;
  DistortionModel * _distortionModel;
  int _LastNTrip; // last triplet count, used to speed up allocation
  double _referenceColor, _sigCol; // average and r.m.s color
  unsigned _nRefrac;
  std::vector<double> _refracCoefficient; // fit parameter
  unsigned int _refracPosInMatrix; // where it stands
  double _JDRef; // average Julian date

  // counts in parameter subsets.
  unsigned int _nParDistortions;
  // unsigned int _nParStars; unused presently
  unsigned int _nParTot;
  unsigned _nMeasuredStars;
  double _posError;  // constant term on error on position (in pixel unit)
  
 public :

  //! this is the only constructor
  AstromFit (Associations &A, DistortionModel *D, double PosError);
  
  //! Does a 1 step minimization, assuming a linear model.
  /*! It calls AssignIndices, LSDerivatives, solves the linear system
    and calls OffsetParams. No line search. Relies on sparse linear
    algebra. */
  unsigned Minimize(const std::string &WhatToFit, const double NSigRejCut=0);

  //! Compute derivatives of measurement terms for this CcdImage
  void LSDerivatives1(const CcdImage &Ccd, 
		      TripletList &TList, Eigen::VectorXd &Rhs,
		      const MeasuredStarList *M=NULL) const;

  //! Compute derivatives of reference terms (if any), associated to the FittedStarList
  void LSDerivatives2(const FittedStarList & Fsl, TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Evaluates the chI^2 derivatives (Jacobian and gradient) for the current WhatToFit setting.
  /*! The Jacobian is provided as triplets, the gradient as a dense
      vector. The parameters which vary are to be set using
      AssignIndices.  */
  void LSDerivatives(TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Set parameter groups fixed or variable and assign indices to each parameter in the big matrix (which will be used by OffsetParams(...).
  void AssignIndices(const std::string &WhatToFit);

  //!The transformations used to propagate errors are freezed to the current state.
  /*! The routine can be called when the mappings are roughly in place.
    After the call, the transformations used to propage errors are no longer
    affected when updating the mappings. This allows to have an exactly linear
    fit, which can be useful. */
  void FreezeErrorScales() {_distortionModel->FreezeErrorScales();}


  //! Offsest the parameters by the requested quantities. The used parameter layout is the one from the last call to AssignIndices or Minimize().
  /*! There is no easy way to check that the current setting of
      WhatToFit and the provided Delta vector are compatible. We can
      only test the size. */
  void OffsetParams(const Eigen::VectorXd &Delta);

  //! Returns a chi2 for the current state
  Chi2 ComputeChi2() const;

  //! Contributions to derivatives from (presumably) outlier terms. No discarding done.
  void OutliersContributions(MeasuredStarList &MOutliers,
			     FittedStarList &FOutLiers,
			     TripletList &TList, 
			     Eigen::VectorXd &Grad);


  //! returns how many outliers were removed. No refit done. MeasOrRef can be "Meas" , "Ref", or "Meas Ref".
  unsigned RemoveOutliers(const double &NSigCut, const std::string &MeasOrRef = "Meas Ref");

  unsigned FindOutliers(const double &NSigCut,
			MeasuredStarList &MSOutliers,
			FittedStarList &FSOutliers,
			const std::string &MeasOrRef = "Meas Ref") const;

  //! Just removes outliers from the fit. No Refit done.
  void RemoveMeasOutliers(MeasuredStarList &Outliers);

  //! Just removes outliers from the fit. No Refit done.
  void RemoveRefOutliers(FittedStarList &Outliers);

  //! Produces both ntuples (cook up names from the provided string)
  void MakeResTuple(const std::string &TupleName) const;

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

  template <class ListType, class Accum> 
    void AccumulateStatImageList(ListType &L, Accum &A) const;

  template <class ImType, class Accum> 
    void AccumulateStatImage(ImType &I, Accum &A) const;

  template <class Accum> 
    void AccumulateStatRefStars(Accum &Accu) const;

  
  //! only for outlier removal
  void GetMeasuredStarIndices(const MeasuredStar &Ms, 
			      std::vector<unsigned> &Indices) const;

};


}}
#endif /* ASTROMFIT__H */
