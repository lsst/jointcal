#include "lsst/meas/simastrom/Projectionhandler.h"
#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/meas/simastrom/CcdImage.h"

namespace lsst {
namespace meas {
namespace simastrom {


class Mapping;


/**********   Stuff for providing Sk22TP gtransfos to a DistortionModel ***/

OneTPPerShoot::OneTPPerShoot(const CcdImageList &L)
{
  for (auto i=L.cbegin(); i!= L.end(); ++i)
    {
      const CcdImage &im = **i;
      if (tMap.find(im.Shoot()) == tMap.end())
	tMap[im.Shoot()] = im.Sky2TP()->Clone();
    }
    
}

const Gtransfo* OneTPPerShoot::Sky2TP(const CcdImage &C) const
{
  auto it=tMap.find(C.Shoot());
  if (it==tMap.end()) return NULL;
  return &*(it->second);
}

}}}
