#ifndef CONSTRAINEDPOLYMODEL__H
#define CONSTRAINEDPOLYMODEL__H

#include "memory" // for std::*_ptr

#include "lsst/jointcal/Eigenstuff.h"

class CcdImage;

#include "lsst/jointcal/DistortionModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
//#include "lsst/jointcal/chiparrangement.h" // just needed for output
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/TwoTransfoMapping.h"
#include "lsst/jointcal/CcdImage.h" // for ShootIdType;


#include <map>


namespace lsst {
namespace jointcal {

class CcdImageList;

//! This is the model used to fit mappings as the combination of a transformation depending on the chip number (instrument model) and a transformation per shoot (anamorphism). The two-transformation Mapping required for this model is TwoTransfoMapping.
/*! This modeling of distortions is meant for set of images from a single mosaic imager. */
class ConstrainedPolyModel : public DistortionModel
{
  /* using ref counts here allows us to not write a destructor nor a copy
     constructor. I could *not* get it to work using std::auto_ptr. */
  typedef std::map<const CcdImage*, std::unique_ptr<TwoTransfoMapping> > mappingMapType;
  mappingMapType _mappings;
  typedef std::map<unsigned, std::unique_ptr<SimpleGtransfoMapping> > chipMapType;
  chipMapType _chipMap;
  typedef std::map<ShootIdType, std::unique_ptr<SimpleGtransfoMapping> > shootMapType;
  shootMapType _shootMap;
  const ProjectionHandler* _sky2TP;
  bool _fittingChips, _fittingShoots;

  Frame _tpFrame; // just for output of the chip transfos
  std::string _instName;

public :
  //!
  ConstrainedPolyModel(const CcdImageList &L,
		       const ProjectionHandler* ProjH,
		       bool InitFromWCS,
		       unsigned NNotFit=0);

  // The following routines are the interface to AstromFit
  //!
  const Mapping* GetMapping(const CcdImage &) const;

  //! Positions the various parameter sets into the parameter vector, starting at FirstIndex
  unsigned AssignIndices(unsigned FirstIndex, std::string &WhatToFit);

  // dispaches the offsets after a fit step into the actual locations of parameters
  void OffsetParams(const Eigen::VectorXd &Delta);

  //! From there on, measurement errors are propagated using the current transfos (and no longer evolve).
  void FreezeErrorScales();

  //! Access to mappings
  const Gtransfo& GetChipTransfo(const unsigned Chip) const;

  //! Access to mappings
  const Gtransfo& GetShootTransfo(const ShootIdType &Shoot) const;

  //! Access to array of shoots involved in the solution.
  std::vector<ShootIdType> GetShoots() const;

  /*! the mapping of sky coordinates (i.e. the coordinate system
  in which fitted stars are reported) onto the Tangent plane
  (into which the pixel coordinates are transformed) */
  const Gtransfo* Sky2TP(const CcdImage &C) const
  { return _sky2TP->Sky2TP(C);}

 //! Cook up a SIP WCS.
  PTR(TanSipPix2RaDec) ProduceSipWcs(const CcdImage &Ccd) const;

  //! Write a transfo file that contains the pixel->tangent plane mappings for each chip.
  /*! These constitute a description of the focal plane
  arrangement. The produced file is used by matchexposure from
  poloka-core */
  //  bool WriteChipArrangement(const std::string &FileName) const;


};

}} // end of namespaces

#endif /* CONSTRAINEDPOLYMODEL__H */
