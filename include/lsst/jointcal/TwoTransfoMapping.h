// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H
#define LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H

#include "memory"

#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimplePolyMapping.h"


namespace lsst {
namespace jointcal {

//! The mapping with two transfos in a row.
class TwoTransfoMapping: public Mapping
{
  SimpleGtransfoMapping *_m1, *_m2;
  unsigned _nPar1, _nPar2;
  struct tmpVars // just there to get around constness issues
  {
    Eigen::MatrixX2d h1,h2;
    Eigen::Matrix2d dt2dx;
  };

  std::unique_ptr<tmpVars> tmp;

  // forbid copies
  TwoTransfoMapping(const TwoTransfoMapping&);
  void operator= (const TwoTransfoMapping&);



 public :
  //!
  TwoTransfoMapping(SimpleGtransfoMapping *ChipM,
		   SimpleGtransfoMapping *VisitM);
  //!
  unsigned Npar() const;

  //!
  void GetMappingIndices(std::vector<unsigned> &Indices) const;

  //!
  void ComputeTransformAndDerivatives(const FatPoint &Where,
				     FatPoint &OutPos,
				     Eigen::MatrixX2d &H)   const;
  //!
  void TransformPosAndErrors(const FatPoint &Where,
			    FatPoint &OutPos) const;

  //!
  void OffsetParams(const double *Delta)
  {// this routine is not used when fitting. used for debugging
    _m1->OffsetParams(Delta);
    _m2->OffsetParams(Delta+ _m1->Npar());
  }

 //! access to transfos
 const Gtransfo& T1() const
 {  return _m1->Transfo();}

 //! access to transfos
 const Gtransfo& T2() const
 {  return _m2->Transfo();}



 //! Currently *not* implemented
 void  PosDerivative(const Point &Where, Eigen::Matrix2d &Der, double  Eps) const;

 //! Currently not implemented
 void FreezeErrorScales();

 private:

 friend class ConstrainedPolyModel;
 //!
 void SetWhatToFit(const bool FittingT1, const bool FittingT2);

};

}} // end of namespaces

#endif // LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H
