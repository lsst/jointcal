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

#ifndef LSST_JOINTCAL_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_ASTROMETRY_MODEL_H

#include <iostream>
#include "memory"

#include "lsst/log/Log.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/AstrometryMapping.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class AstrometryTransform;

/**
 * Interface between AstrometryFit and the combinations of Mappings from pixels to some tangent
 * plane (aka distortions).
 *
 * @param log Logger to send messages to, to keep names consistent when logging.
 */
class AstrometryModel {
public:
    AstrometryModel(LOG_LOGGER log) : _log(log) {}

    /// Return the number of parameters in the mapping of CcdImage
    std::size_t getNpar(CcdImage const &ccdImage) const { return findMapping(ccdImage)->getNpar(); }

    //! Mapping associated to a given CcdImage
    virtual const AstrometryMapping *getMapping(CcdImage const &) const = 0;

    //! Assign indices to parameters involved in mappings, starting at firstIndex. Returns the highest
    //! assigned index.
    virtual Eigen::Index assignIndices(std::string const &whatToFit, Eigen::Index firstIndex) = 0;

    /**
     * Offset the parameters by the provided amounts (by -delta).
     *
     * The shifts are applied according to the indices given in assignIndices.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    //! The transformation used to project the positions of FittedStars.
    /*! This defines the coordinate system into which the Mapping of
        this Ccdimage maps the pixel coordinates. */
    virtual const std::shared_ptr<AstrometryTransform const> getSkyToTangentPlane(
            CcdImage const &ccdImage) const = 0;

    /**
     * Make a SkyWcs that contains this model.
     *
     * @param      ccdImage  The exposure to create the SkyWcs for.
     *
     * @return     SkyWcs containing this model.
     */
    virtual std::shared_ptr<afw::geom::SkyWcs> makeSkyWcs(CcdImage const &ccdImage) const = 0;

    //!
    virtual void freezeErrorTransform() = 0;

    /// Return the total number of parameters in this model.
    virtual std::size_t getTotalParameters() const = 0;

    /**
     * Print a string representation of the contents of this mapping, for debugging.
     *
     * This string representation can be very verbose, as it contains all of the parameters
     * of all of the transforms in this model.
     */
    virtual void print(std::ostream &out) const = 0;

    virtual ~AstrometryModel()= default;;

    /**
     * Return true if this is a "reasonable" model.
     *
     * @param ccdImageList The ccdImages to test the model validity on.
     * @param ndof The number of degrees of freedom in the fit, e.g. from Fitterbase.computeChi2().
     *
     * @return True if the model is valid on all ccdImages.
     */
    bool validate(CcdImageList const &ccdImageList, int ndof) const;

protected:
    /// lsst.logging instance, to be created by a subclass so that messages have consistent name.
    LOG_LOGGER _log;

    /// Return a pointer to the mapping associated with this ccdImage.
    virtual AstrometryMapping *findMapping(CcdImage const &ccdImage) const = 0;
};

std::ostream &operator<<(std::ostream &stream, AstrometryModel const &model);

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTROMETRY_MODEL_H
