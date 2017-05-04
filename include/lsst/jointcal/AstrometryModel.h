// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_ASTROMETRY_MODEL_H

#include "memory"

#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/Mapping.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class Gtransfo;

//! Interface class between AstrometryFit and an actual model for the Mapping (s) from pixels to some tangent plane (aka distortions).
/* For an implementation example, see SimplePolyModel, and the comments at
the top of simplepolymodel.h */
class AstrometryModel
{
public :

  //! Mapping associated to a given CcdImage
  virtual const Mapping* getMapping(const CcdImage &) const = 0;

  //! Assign indices to parameters involved in mappings, starting at FirstIndex. Returns the highest assigned index.
  virtual unsigned assignIndices(unsigned firstIndex, const std::string &whatToFit) = 0;

  //! Offset the parameters by the provided amounts.
  /*! The shifts are applied according to the indices given in
      AssignIndices. */
  virtual void offsetParams(const Eigen::VectorXd &delta) = 0;

  //! The transformation used to project the positions of FittedStars.
  /*! This defines the coordinate system into which the Mapping of
      this Ccdimage maps the pixel coordinates. */
  virtual const Gtransfo* sky2TP(const CcdImage &ccdImage) const = 0;

  //! Cook up a SIP WCS.
  virtual std::shared_ptr<TanSipPix2RaDec> produceSipWcs(const CcdImage &ccdImage) const = 0;

  //!
  virtual void freezeErrorScales() = 0;

  virtual ~AstrometryModel() {};

};


}}

#endif // LSST_JOINTCAL_ASTROMETRY_MODEL_H
