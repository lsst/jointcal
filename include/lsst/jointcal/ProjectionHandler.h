// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PROJECTION_HANDLER_H
#define LSST_JOINTCAL_PROJECTION_HANDLER_H

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Gtransfo.h"
#include "map"

namespace lsst {
namespace jointcal {

class Mapping;
class CcdImage;

/**
 * This is a virtual class that allows a lot of freedom in the choice of the
 * projection from "Sky" (where coodinates are reported) to tangent plane
 * (where they are compared to transformed measurements)
 */
struct ProjectionHandler
{
  virtual const Gtransfo* Sky2TP(const CcdImage &ccdImage) const = 0;

  virtual ~ProjectionHandler() {};

};

/**
 * The simplest implementation of ProjectionHandler. Means that coordinates of
 * objects are expressed in the same space as the arrival mapping space. This
 * is useful for fitting transfo rms between images.
 */
class IdentityProjectionHandler : public ProjectionHandler
{
  GtransfoIdentity id;
 public:
  const Gtransfo* Sky2TP(const CcdImage &ccdImage) const
  { return &id;};

};

/**
 * A projection handler in which all CCDs from the same visit have the same
 * tangent point.
 *
 * We arbitrarily chose that all chips from a same visit have the same tangent
 * point.
 */
class OneTPPerVisitHandler : public ProjectionHandler
{
  typedef std::map<const unsigned, std::shared_ptr<const Gtransfo>> TransfoMap;
  TransfoMap tMap;

 public:
  OneTPPerVisitHandler(const CcdImageList &ccdImageList);

  const Gtransfo* Sky2TP(const CcdImage &ccdImage) const;
};

}} // end of namespaces

#endif // LSST_JOINTCAL_PROJECTION_HANDLER_H
