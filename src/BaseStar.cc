#include <cassert>
#include <fstream>
#include <cstring>  // strstr

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/pex/exceptions.h"


namespace lsst {
namespace jointcal {

bool decreasingFlux(BaseStar const *star1, BaseStar const *star2) {
    return (star1->getFlux() > star2->getFlux());
}

bool increasingMag(BaseStar const *star1, BaseStar const *star2) {
    return (star1->getFlux() < star2->getFlux());
}

/**************** BaseStarList ******************/

int decodeFormat(char const *formatLine, char const *starName) {
    if (!formatLine || !starName) return 0;
    const char *p = strstr(formatLine, starName);
    if (!p) return 0;
    return atoi(p + strlen(starName));
}
}  // namespace jointcal
}  // namespace lsst
