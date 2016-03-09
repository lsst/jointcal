#ifndef PROJECTIONHANDLER__H
#define PROJECTIONHANDLER__H

#include "lsst/jointcal/Gtransfo.h"
#include "map"

namespace lsst {
namespace jointcal {

class Mapping;
class CcdImage;
class CcdImageList;

//! This is a virtual class that allows a lot of freedom in the choice of the projection from "Sky" (where coodinates are reported) to tangent plane (where they are compared to transformed measurements)


struct ProjectionHandler
{
  virtual const Gtransfo* Sky2TP(const CcdImage &C) const = 0;

  virtual ~ProjectionHandler() {};

};

//! the simplest implementation of ProjectionHandler. Means that coordinates of objects are expressed in the same space as the arrival mapping space. This is useful for fitting transfo rms between images.
class IdentityProjection : public ProjectionHandler
{
  GtransfoIdentity id;
 public:
  const Gtransfo* Sky2TP(const CcdImage &C) const
  { return &id;};

};

//! A possible implementation of ProjectionHandler for fitting actual WCS's.
/* ! We arbitrarily chose that all chips from a same shoot have the same tangent point, but other choices can be made. */
class OneTPPerShoot : public ProjectionHandler
{
  typedef std::map<const unsigned, CountedRef<const Gtransfo> > TransfoMap;
  TransfoMap tMap;

 public:
  OneTPPerShoot(const CcdImageList &L);

  const Gtransfo* Sky2TP(const CcdImage &C) const;


};

}} // end of namespaces





#endif
