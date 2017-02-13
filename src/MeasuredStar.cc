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

/******* MeasuredStarList *********/

void  MeasuredStarList::setCcdImage(const CcdImage *ccdImage)
{
  for (auto &i: *this)
      i->setCcdImage(ccdImage);
}

template class StarList<MeasuredStar>; // to force instanciation

}} // end of namespaces






