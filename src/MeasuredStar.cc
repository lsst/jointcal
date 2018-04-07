// -*- C++ -*-
#include <cmath>
#include <vector>

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/jointcal/CcdImage.h"

//#include "preferences.h"
//#include "ccdimage.h"
#include "assert.h"  // for assert

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

BaseStarList &measured2Base(MeasuredStarList &starList) { return (BaseStarList &)starList; }

BaseStarList *measured2Base(MeasuredStarList *starList) { return (BaseStarList *)starList; }

const BaseStarList &measured2Base(const MeasuredStarList &starList) { return (const BaseStarList &)starList; }

const BaseStarList *measured2Base(const MeasuredStarList *starList) { return (BaseStarList *)starList; }

/******* MeasuredStarList *********/

void MeasuredStarList::setCcdImage(const CcdImage *ccdImage) {
    for (auto &i : *this) i->setCcdImage(ccdImage);
}
}  // namespace jointcal
}  // namespace lsst
