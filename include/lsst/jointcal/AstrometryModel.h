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

#include "memory"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/AstrometryMapping.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class Gtransfo;

//! Interface class between AstrometryFit and an actual model for the Mapping (s) from pixels to some tangent
//! plane (aka distortions).
/* For an implementation example, see SimplePolyModel, and the comments at
the top of simplepolymodel.h */
class AstrometryModel {
public:
    /// Return the number of parameters in the mapping of CcdImage
    int getNpar(CcdImage const &ccdImage) const { return findMapping(ccdImage)->getNpar(); }

    //! Mapping associated to a given CcdImage
    virtual const AstrometryMapping *getMapping(CcdImage const &) const = 0;

    //! Assign indices to parameters involved in mappings, starting at firstIndex. Returns the highest
    //! assigned index.
    virtual unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) = 0;

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
    virtual const std::shared_ptr<Gtransfo const> getSky2TP(CcdImage const &ccdImage) const = 0;

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
    virtual int getTotalParameters() const = 0;

    virtual ~AstrometryModel(){};

    /**
     * Return true if this is a "reasonable" model.
     *
     * Not yet implemented: DM-16324
     * An example might be that the model produces finite RA/Dec on each sensor's bounding box.
     *
     * @param ccdImageList The ccdImages to test the model validity on.
     * @return True if the model is valid on all ccdImages.
     */
    bool validate(CcdImageList const &ccdImageList) const { return true; }

protected:
    /// Return a pointer to the mapping associated with this ccdImage.
    virtual AstrometryMapping *findMapping(CcdImage const &ccdImage) const = 0;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTROMETRY_MODEL_H
