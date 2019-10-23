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

#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/CcdImage.h"

namespace lsst {
namespace jointcal {

class Mapping;

std::ostream &operator<<(std::ostream &stream, ProjectionHandler const &projectionHandler) {
    projectionHandler.print(stream);
    return stream;
}

OneTPPerVisitHandler::OneTPPerVisitHandler(const CcdImageList &ccdImageList) {
    for (auto const &i : ccdImageList) {
        const CcdImage &im = *i;
        if (tMap.find(im.getVisit()) == tMap.end()) tMap[im.getVisit()] = im.getSkyToTangentPlane()->clone();
    }
}

const std::shared_ptr<const AstrometryTransform> OneTPPerVisitHandler::getSkyToTangentPlane(
        const CcdImage &ccdImage) const {
    auto it = tMap.find(ccdImage.getVisit());
    if (it == tMap.end()) return nullptr;
    return it->second;
}

void OneTPPerVisitHandler::print(std::ostream &out) const {
    out << "Sky->Tangent Plane projection per visit:" << std::endl;
    for (auto &i : tMap) {
        out << "Visit: " << i.first << std::endl;
        out << *(i.second) << std::endl;
    }
}

}  // namespace jointcal
}  // namespace lsst
