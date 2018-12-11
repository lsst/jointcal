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

#ifndef LSST_JOINTCAL_PROJECTION_HANDLER_H
#define LSST_JOINTCAL_PROJECTION_HANDLER_H

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "map"

namespace lsst {
namespace jointcal {

class Mapping;
class CcdImage;

/**
 * This is a virtual class that allows a lot of freedom in the choice of the
 * projection from "Sky" (where coodinates are reported) to tangent plane
 * (where they are compared to transformed measurements)
 */
struct ProjectionHandler {
    virtual const std::shared_ptr<const AstrometryTransform> getSkyToTangentPlane(
            const CcdImage &ccdImage) const = 0;

    virtual ~ProjectionHandler(){};
};

/**
 * The simplest implementation of ProjectionHandler. Means that coordinates of
 * objects are expressed in the same space as the arrival mapping space. This
 * is useful for fitting transform rms between images.
 */
class IdentityProjectionHandler : public ProjectionHandler {
    std::shared_ptr<AstrometryTransformIdentity> id;

public:
    const std::shared_ptr<const AstrometryTransform> getSkyToTangentPlane(const CcdImage &ccdImage) const {
        return id;
    };
};

/**
 * A projection handler in which all CCDs from the same visit have the same
 * tangent point.
 *
 * We arbitrarily chose that all chips from a same visit have the same tangent
 * point.
 */
class OneTPPerVisitHandler : public ProjectionHandler {
    typedef std::map<const VisitIdType, std::shared_ptr<const AstrometryTransform>> TransformMap;
    TransformMap tMap;

public:
    OneTPPerVisitHandler(const CcdImageList &ccdImageList);

    const std::shared_ptr<const AstrometryTransform> getSkyToTangentPlane(const CcdImage &ccdImage) const;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PROJECTION_HANDLER_H
