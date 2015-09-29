#ifndef DISTORTIONSMODEL__H
#define DISTORTIONSMODEL__H

#include "lsst/meas/simastrom/Eigenstuff.h"
#include "lsst/meas/simastrom/Mapping.h"

namespace lsst {
namespace meas {
namespace simastrom {

class CcdImage;
class Gtransfo;

//! Interface class between AstromFit and an actual model for the Mapping (s) from pixels to some tangent plane (aka distortions).
/* For an implementation example, see SimplePolyModel, and the comments at
the top of simplepolymodel.h */
class DistortionModel
{
public :

  //! Mapping associated to a given CcdImage
  virtual const Mapping* GetMapping(const CcdImage &) const = 0;

  //! Assign indices to parameters involved in mappings, starting at FirstIndex. Returns the highest assigned index.
  virtual unsigned AssignIndices(unsigned FirstIndex, std::string &WhatToFit) = 0;

  //! Offset the parameters by the provided amounts. 
  /*! The shifts are applied according to the indices given in
      AssignIndices. */
  virtual void OffsetParams(const Eigen::VectorXd &Delta) = 0;

  //! The transformation used to project the positions of FittedStars.
  /*! This defines the coordinate system into which the Mapping of
      this Ccdimage maps the pixel coordinates. */
  virtual const Gtransfo* Sky2TP(const CcdImage &C) const = 0;


  //! 
  virtual void FreezeErrorScales() = 0;

};


}}}

#endif /*DISTORTIONSMODEL__H */ 
