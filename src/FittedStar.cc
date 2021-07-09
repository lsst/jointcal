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

#include <iostream>
#include <iomanip>

#include "lsst/log/Log.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/StarList.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.FastFinder");
}

namespace lsst {
namespace jointcal {

// cannot be in fittedstar.h, because of "crossed includes"
FittedStar::FittedStar(const MeasuredStar &measuredStar)
        : BaseStar(measuredStar), _indexInMatrix(-1), _measurementCount(0), _refStar(nullptr) {}

void FittedStar::setRefStar(const RefStar *refStar) {
    if ((_refStar != nullptr) && (refStar != nullptr)) {
        // TODO: should we raise an Exception in this case?
        LOGLS_ERROR(_log, "FittedStar: " << *this
                                         << " is already matched to another RefStar. Clean up your lists.");
        LOGLS_ERROR(_log, "old refStar: " << *_refStar);
        LOGLS_ERROR(_log, "new refStar: " << *refStar);
    } else
        _refStar = refStar;
}

void FittedStar::addMagMeasurement(double magValue, double magWeight) {
    _mag = (_mag * _magErr + magValue * magWeight) / (_magErr + magWeight);
    _magErr += magWeight;
}

/************* FittedStarList ************************/

BaseStarList &Fitted2Base(FittedStarList &This) { return (BaseStarList &)This; }

BaseStarList *Fitted2Base(FittedStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Fitted2Base(const FittedStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Fitted2Base(const FittedStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
