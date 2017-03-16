#ifndef SIMPLEPOLYMODEL__H
#define SIMPLEPOLYMODEL__H

#include "memory"

#include "lsst/jointcal/Eigenstuff.h"

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include <map>

namespace lsst {
namespace jointcal {

class CcdImage;
class CcdImageList;

/* We deal here with coordinate transforms which are fitted
   and/or necessary to AstrometryFit. The classes SimplePolyModel
and SimplePolyMapping implement a model where  there is one
separate transfrom per CcdImage. One could chose other setups.

*/

//! this is the model used to fit independent CCDs, meaning that there is no instrument model.
/* This modeling of distortions can even accommodate images set mixing instruments */
class SimplePolyModel : public AstrometryModel
{
  /* using ref counts here allows us to not write a destructor nor a copy
     constructor. I could *not* get it to work using std::auto_ptr. */
  typedef std::map<const CcdImage*, std::unique_ptr<SimpleGtransfoMapping> > mapType;
  mapType _myMap;
  const ProjectionHandler* _sky2TP;

public :

  //! Sky2TP is just a name, it can be anything
  SimplePolyModel(const CcdImageList &ccdImageList,
		  const ProjectionHandler* projectionHandler,
		  bool initFromWCS,
		  unsigned nNotFit=0,
          unsigned degree=3);

  // The following routines are the interface to AstrometryFit
  //!
  const Mapping* getMapping(const CcdImage &) const;

  //! Positions the various parameter sets into the parameter vector, starting at FirstIndex
  unsigned assignIndices(unsigned firstIndex, const std::string &whatToFit);

  // dispaches the offsets after a fit step into the actual locations of parameters
  void offsetParams(const Eigen::VectorXd &delta);

  /*! the mapping of sky coordinates (i.e. the coordinate system
  in which fitted stars are reported) onto the Tangent plane
  (into which the pixel coordinates are transformed) */
  const Gtransfo* sky2TP(const CcdImage &ccdImage) const
  { return _sky2TP->Sky2TP(ccdImage);}

  //!
  virtual void freezeErrorScales();

  //! Access to mappings
  const Gtransfo& GetTransfo(const CcdImage &ccdImage) const;

  std::shared_ptr<TanSipPix2RaDec> produceSipWcs(const CcdImage &ccdImage) const;

  ~SimplePolyModel() {};
};

}} // end of namespaces

#endif /* SIMPLEPOLYMODEL__H */
