// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <vector>

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/jointcal/CcdImage.h"

//#include "preferences.h"
//#include "ccdimage.h"
#include <cassert>  // for assert

namespace lsst {
namespace jointcal {

/* Interesting fields of the stack catalogs :
'base_SdssCentroid_x'
'base_SdssCentroid_y'
'base_SdssCentroid_xErr'
'base_SdssCentroid_yErr'

We miss the xy uncertainty term.
We can cook it up from the sdss shape:
'base_SdssShape_xx'
'base_SdssShape_yy'
'base_SdssShape_xy'

for fluxes, we might use :
'base_CircularApertureFlux_2_instFlux'
'base_CircularApertureFlux_2_instFluxErr'

 where the '2' should be read from the environment.
*/

BaseStarList &Measured2Base(MeasuredStarList &This) { return (BaseStarList &)This; }

BaseStarList *Measured2Base(MeasuredStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Measured2Base(const MeasuredStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Measured2Base(const MeasuredStarList *This) { return (BaseStarList *)This; }

/******* MeasuredStarList *********/

void MeasuredStarList::setCcdImage(const CcdImage *ccdImage) {
    for (auto &i : *this) i->setCcdImage(ccdImage);
}
}  // namespace jointcal
}  // namespace lsst
