#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/CcdImage.h"

namespace lsst {
namespace jointcal {

class Mapping;

/**********   Stuff for providing Sk22TP gtransfos to a AstrometryModel ***/

OneTPPerVisitHandler::OneTPPerVisitHandler(CcdImageList const &ccdImageList) {
    for (auto const &i : ccdImageList) {
        CcdImage const &im = *i;
        if (tMap.find(im.getVisit()) == tMap.end()) tMap[im.getVisit()] = im.getSky2TP()->clone();
    }
}

Gtransfo const *OneTPPerVisitHandler::getSky2TP(CcdImage const &ccdImage) const {
    auto it = tMap.find(ccdImage.getVisit());
    if (it == tMap.end()) return nullptr;
    return &*(it->second);
}
}  // namespace jointcal
}  // namespace lsst
