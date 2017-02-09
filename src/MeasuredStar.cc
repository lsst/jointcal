// -*- C++ -*-
#include <cmath>
#include <vector>

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/StarList.cc"
#include "lsst/jointcal/CcdImage.h"

//#include "preferences.h"
//#include "ccdimage.h"
#include "assert.h" // for assert


namespace lsst {
namespace jointcal {

  /* Interesting fields of the stack catalogs :
 'base_SdssCentroid_x'
 'base_SdssCentroid_y'
 'base_SdssCentroid_xSigma'
 'base_SdssCentroid_ySigma'

  We miss the xy uncertainty term.
  We can cook it up from the sdss shape:
  'base_SdssShape_xx'
  'base_SdssShape_yy'
  'base_SdssShape_xy'

  for fluxes, we might use :
  'base_CircularApertureFlux_2_flux'
  'base_CircularApertureFlux_2_fluxSigma'

   where the '2' should be read from the environment.
  */

//#define WRITE_STARS

BaseStarList& Measured2Base(MeasuredStarList &This)
{
  return (BaseStarList&) This;
}

BaseStarList* Measured2Base(MeasuredStarList *This)
{
  return (BaseStarList*) This;
}

const BaseStarList& Measured2Base(const MeasuredStarList &This)
{
  return (const BaseStarList &) This;
}

const BaseStarList* Measured2Base(const MeasuredStarList *This)
{
  return (BaseStarList*) This;
}

//! StarList ascii IO's

  std::string MeasuredStar::WriteHeader_(std::ostream & pr ,
					 const char* i) const
{
  //  string format = BaseStar::WriteHeader_(pr,i);
  std::string format = BaseStar::WriteHeader_(pr,i);
  if(i==nullptr) i = "";

  pr << "# eflux" << i << " : " << std::endl
     << "# mag" << i << " : " << std::endl
     << "# wmag" << i << " : " << std::endl
     << "# ccd"<<i<<" : chip number" << std::endl
     << "# expo"<<i<<" : exposure number" << std::endl
     << "# airmass"<<i<<" :  "<<std::endl
     << "# filter"<<i<< " : ugriz=01234 " << std::endl
     << "# valid" << i << " : 0=outlier 1=ok ?" << std::endl;

  format += " MeasuredStar 1 ";
  //  if(i && InitialStarRef)
  //    {
  //      string i2 = "1"+string(i);
  //      format += InitialStarRef->WriteHeader_(pr,i2.c_str());
  //    }
  return format;
}


BaseStar* MeasuredStar::read(std::istream& s, const char* format)
{
  MeasuredStar* ret = new MeasuredStar();
  //  ret->BaseStar::read_it(s, format);
  ret->BaseStar::read_it(s, format);
  double tmp;
  s >> ret->eflux
    >> ret->mag
    >> ret->wmag
    >> tmp // chip
    >> tmp // visit
    >> tmp // airmass
    >> tmp // band
    >> ret->valid;
  return ret;
}

void MeasuredStar::writen(std::ostream& s) const
{
  BaseStar::writen(s);

  s << eflux << ' '
    << (std::isnan(mag)? -1 : mag) << ' '
    << wmag << ' ';

  if(ccdImage)
    std::cout << "ccdImage ! " << ccdImage << std::endl;

  if(ccdImage)
    s << ccdImage->getCcdId() << ' '
      << ccdImage->getVisit() << ' '
      << ccdImage->getAirMass() << ' '
      << ccdImage->getFilter() << ' ';
  else
    s << 0 << ' '
      << 0 << ' '
      << 0 << ' '
      << -1 << ' ';

  s << valid << ' ';

}

/******* MeasuredStarList *********/

void  MeasuredStarList::setCcdImage(const CcdImage *C)
{
  for (MeasuredStarIterator i= begin(); i != end(); ++i)
      (*i)->setCcdImage(C);
}

template class StarList<MeasuredStar>; // to force instanciation

}} // end of namespaces






