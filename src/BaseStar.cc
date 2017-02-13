#include <assert.h>
#include <fstream>
#include <string.h> // strstr

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.cc"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

static double sq(double x) { return x*x;}

bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2)
{
return (S1->flux > S2->flux);
}

bool IncreasingMag(const BaseStar *S1, const BaseStar *S2)
{
return (S1->flux < S2->flux);
}


/**************** BaseStarList ******************/

int DecodeFormat(const char *FormatLine, const char *StarName)
{
if (!FormatLine || !StarName) return 0;
const char *p= strstr(FormatLine, StarName);
 if (!p) return  0;
return atoi( p + strlen(StarName));
}


  //Instanciate all routines from the template
template class StarList<BaseStar>; /* to force instanciation */


}} // end of namespaces
