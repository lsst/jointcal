#include <assert.h>
#include <fstream>
#include <string.h>  // strstr

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

bool decreasingFlux(const BaseStar *star1, const BaseStar *star2) {
    return (star1->getFlux() > star2->getFlux());
}

bool increasingMag(const BaseStar *star1, const BaseStar *star2) {
    return (star1->getFlux() < star2->getFlux());
}

/**************** BaseStarList ******************/

int decodeFormat(const char *formatLine, const char *starName) {
    if (!formatLine || !starName) return 0;
    const char *p = strstr(formatLine, starName);
    if (!p) return 0;
    return atoi(p + strlen(starName));
}
}  // namespace jointcal
}  // namespace lsst
