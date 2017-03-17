#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/CcdImage.h"

namespace lsst {
namespace jointcal {


class Mapping;


/**********   Stuff for providing Sk22TP gtransfos to a AstrometryModel ***/

OneTPPerVisitHandler::OneTPPerVisitHandler(const CcdImageList &ccdImageList)
{
  for (auto const &i: ccdImageList)
    {
      const CcdImage &im = *i;
      if (tMap.find(im.getVisit()) == tMap.end())
	tMap[im.getVisit()] = im.Sky2TP()->Clone();
    }

}

const Gtransfo* OneTPPerVisitHandler::Sky2TP(const CcdImage &ccdImage) const
{
  auto it=tMap.find(ccdImage.getVisit());
  if (it==tMap.end()) return nullptr;
  return &*(it->second);
}

}} // end of namespaces
