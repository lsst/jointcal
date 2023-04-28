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
#include <iostream>
#include "lsst/geom/SpherePoint.h"
#include "lsst/jointcal/ProperMotion.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

Point ProperMotion::apply(const Point& star, double timeDeltaYears) const {
    geom::SpherePoint spherePoint(star.x, star.y, geom::degrees);
    double amount = std::hypot(_ra * timeDeltaYears, _dec * timeDeltaYears);
    // If delta-time is negative, the correction is in the opposite direction.
    amount = timeDeltaYears < 0 ? -amount : amount;
    auto result = spherePoint.offset(_offsetBearing * geom::radians, amount * geom::radians);
    Point newStar(star);
    newStar.x = result.getRa().asDegrees();
    newStar.y = result.getDec().asDegrees();
    return newStar;
}

std::ostream &operator<<(std::ostream &stream, ProperMotion const &pm) {
    stream << "pm_ra*cos(dec)=" << pm._ra << "rad/yr, pm_dec=" << pm._dec << "rad/yr, pm_raErr=" << pm._raErr
           << "rad/yr, pm_decErr=" << pm._decErr << "rad/yr, pm_raDecCov=" << pm._raDecCov;
    return stream;
}

}  // namespace jointcal
}  // namespace lsst
