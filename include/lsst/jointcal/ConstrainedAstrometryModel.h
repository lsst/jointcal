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

#ifndef LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H

#include "memory"  // for std::*_ptr

#include "lsst/jointcal/Eigenstuff.h"

class CcdImage;

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/ChipVisitAstrometryMapping.h"
#include "lsst/jointcal/CcdImage.h"

#include <map>

namespace lsst {
namespace jointcal {

/**
 * A multi-component model, fitting mappings for sensors and visits simultaneously.
 *
 * This is the model used to fit mappings as the combination of a
 * transformation depending on the chip number (instrument model) and a
 * transformation per visit (anamorphism). The two-transformation Mapping
 * required for this model is ChipVisitAstrometryMapping. This modeling of distortions
 * is meant for a set of images from a single mosaic imager.
 *
 * @param ccdImageList The exposures that will be fit.
 * @param projectionHandler The projection from "Sky" (where the "true" coordinates live) to "Tangent Plane"
 *                          (where the fitting occurs).
 * @param chipOrder The polynomial order of the pixel->focal plane mapping for each sensor.
 * @param visitOrder The polynomial order of the focal plane->tangent plane mapping for each visit.
 */
class ConstrainedAstrometryModel : public AstrometryModel {
public:
    ConstrainedAstrometryModel(CcdImageList const &ccdImageList,
                               std::shared_ptr<ProjectionHandler const> projectionHandler, int chipOrder,
                               int visitOrder);

    /// No copy or move: there is only ever one instance of a given model (i.e. per ccd+visit)
    ConstrainedAstrometryModel(ConstrainedAstrometryModel const &) = delete;
    ConstrainedAstrometryModel(ConstrainedAstrometryModel &&) = delete;
    ConstrainedAstrometryModel &operator=(ConstrainedAstrometryModel const &) = delete;
    ConstrainedAstrometryModel &operator=(ConstrainedAstrometryModel &&) = delete;

    // The following routines are the interface to AstrometryFit
    //!
    AstrometryMapping const *getMapping(CcdImage const &) const override;

    /**
     * Positions the various parameter sets into the parameter vector, starting at
     * firstIndex.
     */
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) override;

    /**
     * Dispaches the offsets after a fit step into the actual locations of
     * parameters.
     */
    void offsetParams(Eigen::VectorXd const &Delta) override;

    /**
     * From there on, measurement errors are propagated using the current
     * transforms (and no longer evolve).
     */
    void freezeErrorTransform() override;

    /// @copydoc AstrometryModel::getTotalParameters
    int getTotalParameters() const override;

    //! Access to mappings
    AstrometryTransform const &getChipTransform(CcdIdType const chip) const;

    //! Access to mappings
    AstrometryTransform const &getVisitTransform(VisitIdType const &visit) const;

    //! Access to array of visits involved in the solution.
    std::vector<VisitIdType> getVisits() const;

    /**
     * The mapping of sky coordinates (i.e. the coordinate system in which fitted
     * stars are reported) onto the Tangent plane (into which the pixel coordinates
     * are transformed).
     */
    const std::shared_ptr<AstrometryTransform const> getSkyToTangentPlane(
            CcdImage const &ccdImage) const override {
        return _skyToTangentPlane->getSkyToTangentPlane(ccdImage);
    }

    /// @copydoc AstrometryModel::makeSkyWcs
    std::shared_ptr<afw::geom::SkyWcs> makeSkyWcs(CcdImage const &ccdImage) const override;

private:
    std::unordered_map<CcdImageKey, std::unique_ptr<ChipVisitAstrometryMapping>> _mappings;
    std::map<CcdIdType, std::shared_ptr<SimpleAstrometryMapping>> _chipMap;
    std::map<VisitIdType, std::shared_ptr<SimpleAstrometryMapping>> _visitMap;
    const std::shared_ptr<ProjectionHandler const> _skyToTangentPlane;
    bool _fittingChips, _fittingVisits;

    /// @copydoc AstrometryModel::findMapping
    AstrometryMapping *findMapping(CcdImage const &ccdImage) const override;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H
