#ifndef PHOTOMMODEL__H
#define PHOTOMMODEL__H

#include "lsst/meas/simastrom/Eigenstuff.h"


namespace lsst {
namespace meas {
namespace simastrom {

class CcdImage;
class Point;
class MeasuredStar;

//! Interface class for PhotomFit
class PhotomModel
{
public :

  //! Assign indices to parameters involved in mappings, starting at FirstIndex. Returns the highest assigned index.
  virtual unsigned AssignIndices(unsigned FirstIndex) = 0;

  //! Offset the parameters by the provided amounts. 
  /*! The shifts are applied according to the indices given in
      AssignIndices. */
  virtual void OffsetParams(const Eigen::VectorXd &Delta) = 0;

  virtual double PhotomFactor(const Point &Where, const CcdImage& C) const =0;

  virtual unsigned GetIndicesAndDerivatives(const MeasuredStar &M,
					    const CcdImage &Ccd, 
					      std::vector<unsigned> &Indices,
					      Eigen::VectorXd &D) = 0;


  virtual ~PhotomModel() {};

};


}}}

#endif /*DISTORTIONSMODEL__H */ 
