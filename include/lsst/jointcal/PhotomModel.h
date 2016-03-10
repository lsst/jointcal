#ifndef PHOTOMMODEL__H
#define PHOTOMMODEL__H

#include "lsst/jointcal/Eigenstuff.h"
#include <string>
#include <vector>

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;
class MeasuredStar;

//! Interface class for PhotomFit
class PhotomModel
{
public :

  //! Assign indices to parameters involved in mappings, starting at FirstIndex. Returns the highest assigned index.
  virtual unsigned AssignIndices(const std::string &WhatToFit, unsigned FirstIndex) = 0;

  //! Offset the parameters by the provided amounts.
  /*! The shifts are applied according to the indices given in
      AssignIndices. */
  virtual void OffsetParams(const Eigen::VectorXd &Delta) = 0;

  //! Where is to be expressed in Ccd coordinates.
  virtual double PhotomFactor(const CcdImage& C, const Point &Where) const =0;

  //! number of parameters to be read in Indices.size()
  virtual void GetIndicesAndDerivatives(const MeasuredStar &M,
					const CcdImage &Ccd,
					std::vector<unsigned> &Indices,
					Eigen::VectorXd &D) = 0;


  virtual ~PhotomModel() {};

};


}} // end of namespaces

#endif /*DISTORTIONSMODEL__H */
