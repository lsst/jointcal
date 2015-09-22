#ifndef ASTROMFIT__H
#define ASTROMFIT__H

#include <string>
#include <iostream>
#include <sstream>
#include <map>

#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Eigenstuff.h"
#include "lsst/meas/simastrom/Tripletlist.h"
#include "lsst/meas/simastrom/DistortionModel.h"

namespace lsst {
namespace meas {
namespace simastrom {

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

 struct Chi2;

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
  bool Minimize(const std::string &WhatToFit);

  //! Compute the derivatives of the measurement terms for this CcdImage
  void LSDerivatives1(const CcdImage &Ccd, 
		     TripletList &TList, Eigen::VectorXd &Rhs) const;

  //! Compute the derivatives of the reference terms
  void LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) const;

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

  //! returns how many outliers were removed. No refit done.
  unsigned RemoveOutliers(const double &NSigCut);

  //! Produces a tuple containing residuals of measurement terms.
  void MakeMeasResTuple(const std::string &TupleName) const;

  //! Produces a tuple containing residuals of reference terms.
  void MakeRefResTuple(const std::string &TupleName) const;

  //! Produces both ntuples (cook up names from the provided string)
  void MakeResTuple(const std::string &TupleName) const;

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
  
  //! only for outlier removal
  void GetMeasuredStarIndices(const MeasuredStar &Ms, 
			      std::vector<unsigned> &Indices) const;

};

//! Simple structure to accumulate Chi2 and Ndof
struct Chi2
{
  double chi2;
  unsigned ndof;

  Chi2() : chi2(0), ndof(0) {};

  friend std::ostream& operator << (std::ostream& s, const Chi2 &C)
    {
      s << "Chi2/ndof : " << C.chi2 << '/' << C.ndof << '=' <<  C.chi2/C.ndof; return s;
    }

  std::string __str__()
  {
    std::stringstream s;
    s << "Chi2/ndof : " << chi2 << '/' << ndof << '=' <<  chi2/ndof; 
    return s.str();
  }

  // Addentry has a third argument in order to make it compatible with an 
  //other stat accumulator.
  void AddEntry(double Inc, unsigned Dof, const MeasuredStar *M) 
  {chi2+= Inc; ndof += Dof;}

  void operator += (const Chi2 &R) 
  {chi2 += R.chi2; ndof += R.ndof;}

}; // end of struct Chi2




}}}
#endif /* ASTROMFIT__H */
