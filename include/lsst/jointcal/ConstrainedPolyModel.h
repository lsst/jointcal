#ifndef CONSTRAINEDPOLYMODEL__H
#define CONSTRAINEDPOLYMODEL__H

#include "memory" // for std::*_ptr

#include "lsst/jointcal/Eigenstuff.h"

class CcdImage;

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/TwoTransfoMapping.h"
#include "lsst/jointcal/CcdImage.h" // for VisitIdType;


#include <map>


namespace lsst {
namespace jointcal {

class CcdImageList;

/**
 * This is the model used to fit mappings as the combination of a
 * transformation depending on the chip number (instrument model) and a
 * transformation per visit (anamorphism). The two-transformation Mapping
 * required for this model is TwoTransfoMapping. This modeling of distortions
 * is meant for set of images from a single mosaic imager.
 */
class ConstrainedPolyModel : public AstrometryModel
{
  // NOTE: Using ref counts here allows us to not write a destructor nor a copy
  // constructor. I could *not* get it to work using std::auto_ptr.
  typedef std::map<const CcdImage*, std::unique_ptr<TwoTransfoMapping> > mappingMapType;
  mappingMapType _mappings;
  typedef std::map<unsigned, std::unique_ptr<SimpleGtransfoMapping> > chipMapType;
  chipMapType _chipMap;
  typedef std::map<VisitIdType, std::unique_ptr<SimpleGtransfoMapping> > visitMapType;
  visitMapType _visitMap;
  const ProjectionHandler* _sky2TP;
  bool _fittingChips, _fittingVisits;

  Frame _tpFrame; // just for output of the chip transfos

public :
  ConstrainedPolyModel(const CcdImageList &ccdImageList,
                       const ProjectionHandler* projectionHandler,
                       bool initFromWCS,
                       unsigned nNotFit=0);

  // The following routines are the interface to AstrometryFit
  //!
  const Mapping* getMapping(const CcdImage &) const;

  /**
   * Positions the various parameter sets into the parameter vector, starting at
   * FirstIndex.
   */
  unsigned assignIndices(unsigned FirstIndex, std::string &WhatToFit);

  /**
   * Dispaches the offsets after a fit step into the actual locations of
   * parameters.
   */
  void offsetParams(const Eigen::VectorXd &Delta);

  /**
   * From there on, measurement errors are propagated using the current
   * transfos (and no longer evolve).
   */
  void freezeErrorScales();

  //! Access to mappings
  const Gtransfo& getChipTransfo(const unsigned Chip) const;

  //! Access to mappings
  const Gtransfo& getVisitTransfo(const VisitIdType &Visit) const;

  //! Access to array of visits involved in the solution.
  std::vector<VisitIdType> getVisits() const;

  /**
   * The mapping of sky coordinates (i.e. the coordinate system in which fitted
   * stars are reported) onto the Tangent plane (into which the pixel coordinates
   * are transformed).
   */
  const Gtransfo* sky2TP(const CcdImage &ccdImage) const
  { return _sky2TP->Sky2TP(ccdImage);}

  std::shared_ptr<TanSipPix2RaDec> produceSipWcs(const CcdImage &ccdImage) const;
};

}} // end of namespaces

#endif /* CONSTRAINEDPOLYMODEL__H */
