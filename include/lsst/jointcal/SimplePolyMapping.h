#ifndef SIMPLEPOLYMAPPING__H
#define SIMPLEPOLYMAPPING__H

#include <memory> // for unique_ptr

#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/CcdImage.h"

//! Class for a simple mapping implementing a generic Gtransfo
/*! It uses a template rather than a pointer so that the derived
classes can use the specifics of the transfo. The class
simplePolyMapping overloads a few routines. */

namespace lsst {
namespace jointcal {

class SimpleGtransfoMapping : public Mapping
{

 protected:
  bool toFit;
  unsigned index;
  /* inheritance may also work. Perhaps with some trouble because
     some routines in Mapping and Gtransfo have the same name */
  std::shared_ptr<Gtransfo> transfo;

  std::shared_ptr<Gtransfo> errorProp;
  /* to avoid allocation at every call of PosDerivatives.
     use a pointer for constness */
  std::unique_ptr<GtransfoLin> lin;


#ifdef STORAGE
  //! this is modern compilation-time check:
  static_assert(std::is_base_of<Gtransfo, Tr>::value,
		"SimpleGtransfoMapping requires a class deriving from Gtransfo");
#endif


 public :

 SimpleGtransfoMapping(const Gtransfo &T, bool ToFit=true) : toFit(ToFit), transfo(T.Clone()), errorProp(transfo), lin(new GtransfoLin)
  {
    // in this order:
    // take a copy of the input transfo,
    // assign the transformation used to propagate errors to the transfo itself
    // reserve some memory space to compute the derivatives (efficiency).
  }

  virtual void FreezeErrorScales()
  {
    // from there on, updating the transfo does not change the errors.
    errorProp.reset(transfo->Clone());
  }

  // interface Mapping functions:

  //!
  unsigned Npar() const
  {
    if (toFit) return transfo->Npar();
    else return 0;
  }

  //!
  void GetMappingIndices(std::vector<unsigned> &Indices) const
  {
    if (Indices.size() < Npar()) Indices.resize(Npar());
    for (unsigned k=0; k<Npar(); ++k) Indices[k] = index+k;
  }

  //!
  void TransformPosAndErrors(const FatPoint &Where, FatPoint &OutPos) const
  {
    transfo->TransformPosAndErrors(Where,OutPos);
    FatPoint tmp;
    errorProp->TransformPosAndErrors(Where,tmp);
    OutPos.vx = tmp.vx;
    OutPos.vy = tmp.vy;
    OutPos.vxy = tmp.vxy;
  }

  //!
  void  PosDerivative(const Point &Where, Eigen::Matrix2d &Der, double Eps) const
  {
    errorProp->Derivative(Where, *lin, Eps);
    Der(0,0) = lin->Coeff(1,0,0);
    //
    /* This does not work : it was proved by rotating the frame
       see the compilation switch ROTATE_T2 in constrainedpolymodel.cc
    Der(1,0) = lin->Coeff(1,0,1);
    Der(0,1) = lin->Coeff(0,1,0);
    */
    Der(1,0) = lin->Coeff(0,1,0);
    Der(0,1) = lin->Coeff(1,0,1);
    Der(1,1) = lin->Coeff(0,1,1);
  }

  //!
  void OffsetParams(const double *Delta)
  {
    transfo->OffsetParams(Delta);
  }

  //! position of the parameters within the grand fitting scheme
  unsigned Index() const { return index;}

  //!
  void  SetIndex(unsigned I) {index=I;}

  virtual void ComputeTransformAndDerivatives(const FatPoint &Where,
                                              FatPoint &OutPos,
                                              Eigen::MatrixX2d &H) const
  {
    TransformPosAndErrors(Where,OutPos);
    transfo->ParamDerivatives(Where, &H(0,0), &H(0,1));
  }

  //! Access to the (fitted) transfo
  virtual const Gtransfo&  Transfo() const {return *transfo;}

};

//! Mapping implementation for a polynomial transformation.
class SimplePolyMapping : public SimpleGtransfoMapping
{
  /* to better condition the 2nd derivative matrix, the
  transformed coordinates are mapped (roughly) on [-1,1].
  We need both the transform and its derivative. */
  GtransfoLin _centerAndScale;
  Eigen::Matrix2d preDer;

  /* Where we store the combination. We use a pointer for
  constness. Could not get it to work with smart pointers.
  */
  GtransfoPoly* actualResult;

 public:

  ~SimplePolyMapping() { delete actualResult;}



  // ! contructor.
  /*! The transformation will be initialized to Transfo, so that the effective transformation
    reads Transfo*CenterAndScale */
 SimplePolyMapping(const GtransfoLin& CenterAndScale, const GtransfoPoly& Transfo):
  SimpleGtransfoMapping(Transfo), _centerAndScale(CenterAndScale)
  {
    // We assume that the initialization was done properly, for example that
    // Transfo = pix2TP*CenterAndScale.invert(), so we do not touch transfo.
    /* store the (spatial) derivative of _centerAndScale. For the extra
       diagonal terms, just copied the ones in PosDerivatives */
    preDer(0,0) = _centerAndScale.Coeff(1,0,0);
    preDer(1,0) = _centerAndScale.Coeff(0,1,0);
    preDer(0,1) = _centerAndScale.Coeff(1,0,1);
    preDer(1,1) = _centerAndScale.Coeff(0,1,1);

    // reserve space for the result
    actualResult = new GtransfoPoly();

    // check of matrix indexing (once for all)
    MatrixX2d H(3,2);
    assert((&H(1,0) - &H(0,0)) == 1);
  }

  /* The SimpleGtransfoMapping version does not account for the
     _centerAndScale transfo */

  void  PosDerivative(const Point &Where, Eigen::Matrix2d &Der,
		      double  Eps) const
  {
    Point tmp = _centerAndScale.apply(Where);
    errorProp->Derivative(tmp, *lin, Eps);
    Der(0,0) = lin->Coeff(1,0,0);
    //
    /* This does not work : it was proved by rotating the frame
       see the compilation switch ROTATE_T2 in constrainedpolymodel.cc
    Der(1,0) = lin->Coeff(1,0,1);
    Der(0,1) = lin->Coeff(0,1,0);
    */
    Der(1,0) = lin->Coeff(0,1,0);
    Der(0,1) = lin->Coeff(1,0,1);
    Der(1,1) = lin->Coeff(0,1,1);
    Der = preDer*Der;
  }



  //! Calls the transforms and implements the centering and scaling of coordinates
  /* We should put the computation of error propagation and
     parameter derivatives into the same Gtransfo routine because
     it could be significantly faster */
  virtual void ComputeTransformAndDerivatives(const FatPoint &Where,
                                              FatPoint &OutPos,
                                              Eigen::MatrixX2d &H) const
    {
      FatPoint mid;
      _centerAndScale.TransformPosAndErrors(Where,mid);
      transfo->TransformPosAndErrors(mid,OutPos);
      FatPoint tmp;
      errorProp->TransformPosAndErrors(mid,tmp);
      OutPos.vx = tmp.vx;
      OutPos.vy = tmp.vy;
      OutPos.vxy = tmp.vxy;
      transfo->ParamDerivatives(mid, &H(0,0), &H(0,1));
    }

  //! Implements as well the centering and scaling of coordinates
  void TransformPosAndErrors(const FatPoint &Where, FatPoint &OutPos) const
  {
    FatPoint mid;
    _centerAndScale.TransformPosAndErrors(Where,mid);
    transfo->TransformPosAndErrors(mid,OutPos);
    FatPoint tmp;
    errorProp->TransformPosAndErrors(mid,tmp);
    OutPos.vx = tmp.vx;
    OutPos.vy = tmp.vy;
    OutPos.vxy = tmp.vxy;
  }

  //! Access to the (fitted) transfo
  const Gtransfo& Transfo() const
  {
    // Cannot fail given the contructor:
    const GtransfoPoly *fittedPoly = dynamic_cast<const GtransfoPoly*>(&(*transfo));
    *actualResult = (*fittedPoly)*_centerAndScale;
    return *actualResult;
  }


};


#ifdef STORAGE
/*! "do nothing" mapping. The Ccdimage's that "use" this one impose the
  coordinate system */
class SimpleIdentityMapping : public SimpleGtransfoMapping<GtransfoIdentity>
{
 public:

  //! nothing to do.
  virtual void ComputeTransformAndDerivatives(const FatPoint &Where,
                                              FatPoint &OutPos,
                                              Eigen::MatrixX2d &H) const
  {
    OutPos = Where;
  }

};
#endif

}} // end of namespaces
#endif
