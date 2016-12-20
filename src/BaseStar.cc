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


void BaseStar::read_it(std::istream & rd, const char *format)
{
 int formatValue = 0;
 if (format)
   formatValue = DecodeFormat(format,"BaseStar");
 if (formatValue == 3)
   {
     rd >> x >> y >> vx >> vy >> vxy >> flux;
     /* write sig(x) sig(y), rho ... */
     vxy *= (vx*vy);
     vx *= vx;
     vy *= vy;
   }
 else if (formatValue == 2 || formatValue == 0)
   rd >> x >> y >> flux;
 else if (formatValue == 1) // only shift back if shifted when written
   {
     rd >> x >> y >> flux;
     x -= MEMPIX2DISK;
     y -= MEMPIX2DISK;
   }
 else throw  LSST_EXCEPT(pexExcept::OutOfRangeError," Unknown format value for BaseStar ");
}

/* this routine is used by DicStar. It *HAS* to be consistent with BaseStar::read_it, above */
unsigned NValsBaseStar(const char *Format)
{
  int formatValue = DecodeFormat(Format, "BaseStar");
  if (formatValue <= 2) return 3;
  if (formatValue == 3) return 6;
  return 0;
}



BaseStar* BaseStar::read(std::istream & rd, const char *format)
{
  BaseStar *p = new BaseStar();
  p->read_it(rd, format);
  return p;
}


std::string BaseStar::WriteHeader_(std::ostream & stream, const char*i) const
{
  if (i==nullptr) i = "";
  stream << "# x"<< i <<" : x position (pixels)" << std::endl
	 << "# y"<< i <<" : y position (pixels)" << std::endl
	 << "# sx"<< i <<" : x position r.m.s " << std::endl
	 << "# sy"<< i <<" : y position r.m.s " << std::endl
	 << "# rhoxy"<< i <<" : xy correlation " << std::endl
	 << "# flux"<< i <<" : flux in image ADUs" << std::endl ;
  return " BaseStar 3 ";
}

void BaseStar::WriteHeader(std::ostream & stream) const
{
  std::string format = WriteHeader_(stream);
  stream << "# format " << format << std::endl;
  stream << "# end " << std::endl ;
};


void BaseStar::writen(std::ostream &s) const
{
  assert(vx>0 && vy>0 && sq(vxy)<vx*vy);
  /* write (sigx,sigy,rho) rather than (vx,vy,vxy).
     This limits shortcomings of truncation, and is more useful
     when reading or plotting.
  */
  double sx = sqrt(vx);
  double sy = sqrt(vy);
  s << x << ' ' << y << ' '
    << sx << ' ' << sy << ' ' << vxy/sx/sy << ' '
    << flux << ' ' ;
}


void BaseStar::write(std::ostream &s) const
{
  writen(s);
  s << std::endl ;
}

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
