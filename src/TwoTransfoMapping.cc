#include "lsst/jointcal/TwoTransfoMapping.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;


namespace lsst {
namespace jointcal {


TwoTransfoMapping::TwoTransfoMapping(SimpleGtransfoMapping *M1, 
				     SimpleGtransfoMapping *M2): 
  _m1(M1), _m2(M2)
{
  /* Allocate the record of temporary variables, so that they are not
     allocated at every call. This is hidden behind a pointer in order
     to be allowed to alter them in a const routine. */
  tmp = std::unique_ptr<tmpVars>(new tmpVars);
  SetWhatToFit(true,true);
}

unsigned TwoTransfoMapping::Npar() const
{
  return _nPar1+ _nPar2;
}


void TwoTransfoMapping::GetMappingIndices(std::vector<unsigned> &Indices) const
{
  unsigned npar = Npar();
  if (Indices.size()<npar) Indices.resize(npar); 
  // in case we are only fitting one of the two transfos
  if (_nPar1) _m1->GetMappingIndices(Indices);
  else if (_nPar2) {_m2->GetMappingIndices(Indices); return;}
  // if we get here we are fitting both
  // there is probably a more elegant way to feed a subpart of a std::vector
  std::vector<unsigned> ind2(_nPar2);
  _m2->GetMappingIndices(ind2);
  for (unsigned k=0; k<_nPar2; ++k) Indices.at(k+_nPar1) = ind2.at(k);
}

void TwoTransfoMapping::ComputeTransformAndDerivatives(const FatPoint &Where,
						       FatPoint &OutPos,
						       Eigen::MatrixX2d &H) 
  const
{
  // not true in general. Will crash if H is too small.
  //  assert(H.cols()==Npar());

  FatPoint pMid;
  // don't need errors there but no Mapping::Transform() routine.

  if (_nPar1)
    {
      _m1->ComputeTransformAndDerivatives(Where, pMid,  tmp->h1);
      // the last argument is epsilon and is not used for polynomials      
      _m2->PosDerivative(pMid, tmp->dt2dx, 1e-4); 
      H.block(0,0,_nPar1,2) = tmp->h1*tmp->dt2dx;
    }
  else _m1->TransformPosAndErrors(Where, pMid);
  if (_nPar2)
    {
      _m2->ComputeTransformAndDerivatives(pMid,OutPos, tmp->h2);
      H.block(_nPar1,0,_nPar2,2) = tmp->h2;
    }
  else _m2->TransformPosAndErrors(pMid, OutPos);

}

  /*! Sets the _nPar{1,2} and allocates H matrices accordingly, to
     avoid allocation at every call. If we did not care about dynamic
     allocation, we could just put the information of what moves and
     what doesn't into the SimpleGtransfoMapping. */
void TwoTransfoMapping::SetWhatToFit(const bool FittingT1, const bool FittingT2)
{

  if (FittingT1)
    {
      _nPar1 = _m1->Npar();
      tmp->h1 = Eigen::MatrixX2d(_nPar1,2);
    }
  else _nPar1 = 0;
  if (FittingT2)
    {
      _nPar2 = _m2->Npar();
      tmp->h2 = Eigen::MatrixX2d(_nPar2,2);
    }
  else _nPar2 = 0;
}


void TwoTransfoMapping::TransformPosAndErrors(const FatPoint &Where,
					      FatPoint &OutPos) const
{
  FatPoint pMid;
  _m1->TransformPosAndErrors(Where, pMid);
  _m2->TransformPosAndErrors(pMid, OutPos);
}

void  TwoTransfoMapping::PosDerivative(const Point &Where, 
				       Eigen::Matrix2d &Der, 
				       const double & Eps) const
{
  Eigen::Matrix2d d1,d2; // seems that it does not trigger dynamic allocation
  _m1->PosDerivative(Where, d1, 1e-4);
  FatPoint pMid;
  _m1->TransformPosAndErrors(Where, pMid);
  _m2->PosDerivative(pMid, d2, 1e-4);
  /* The following line is not a mistake. It is a consequence 
     of chosing Der(0,1) = d(y_out)/d x_in. */
  Der = d1*d2;
} 

void TwoTransfoMapping::FreezeErrorScales()
{
  throw LSST_EXCEPT(pexExcept::TypeError," The routine  TwoTransfoMapping::FreezeErrorScales() was thought to be useless and is not implemented (yet)"); 
} 

}} // end of namespaces

