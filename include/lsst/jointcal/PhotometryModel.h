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

#ifndef LSST_JOINTCAL_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_PHOTOMETRY_MODEL_H

#include <string>
#include <vector>

#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/log/Log.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {

class PhotometryModel {
public:
    /**
     * @param log Logger to send messages to, to keep names consistent when logging.
     * @param errorPedestal_ Pedestal on flux/magnitude error (percent of flux or delta magnitude).
     */
    PhotometryModel(LOG_LOGGER log, double errorPedestal_ = 0) : _log(log), errorPedestal(errorPedestal_) {}

    /**
     * Assign indices in the full matrix to the parameters being fit in the mappings, starting at firstIndex.
     *
     * @param[in]  whatToFit   String containing parameters to fit.
     * @param[in]  firstIndex  Index to start assigning at.
     *
     * @return     The highest assigned index.
     */
    virtual unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) = 0;

    /**
     * Offset the parameters by the provided amounts (by -delta).
     *
     * The shifts are applied according to the indices given in assignIndices.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    /**
     * Offset the appropriate flux or magnitude (by -delta).
     *
     * @param fittedStar The star to update.
     * @param delta The amount to update by.
     */
    virtual void offsetFittedStar(FittedStar &fittedStar, double delta) const = 0;

    /**
     * Compute the residual between the model applied to a star and its associated fittedStar.
     *
     * @f[
     *     residual = Model(measuredStar) - fittedStar
     * @f]
     *
     * @param ccdImage The ccdImage where measuredStar resides.
     * @param measuredStar The measured star position to compute the residual of.
     *
     * @return The residual.
     */
    virtual double computeResidual(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const = 0;

    /**
     * Return the on-sky transformed flux for measuredStar on ccdImage.
     *
     * @param[in]  ccdImage     The ccdImage where measuredStar resides.
     * @param      measuredStar The measured star position to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
     */
    virtual double transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const = 0;

    /**
     * Return the on-sky transformed flux uncertainty for measuredStar on ccdImage.
     * Identical to transform() until freezeErrorTransform() is called.
     *
     * @param[in]  ccdImage     The ccdImage where measuredStar resides.
     * @param      measuredStar The measured star position to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
     */
    virtual double transformError(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const = 0;

    /**
     * Once this routine has been called, the error transform is not modified by offsetParams().
     *
     * The routine can be called when the mappings are roughly in place. After the call, the transformations
     * used to propagate errors are no longer affected when updating the mappings. This allows an exactly
     * linear fit, which can be necessary for some model+data combinations.
     */
    virtual void freezeErrorTransform() = 0;

    /**
     * Get how this set of parameters (of length Npar()) map into the "grand" fit.
     *
     * @param[in]  ccdImage  The ccdImage to look up.
     * @param[out] indices   The indices of the mapping associated with ccdImage.
     */
    virtual void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) const = 0;

    /**
     * Compute the parametric derivatives of this model.
     *
     * @param[in]   measuredStar  The measured star with the position and flux to compute at.
     * @param[in]   ccdImage      The ccdImage containing the measured star, to find the correct mapping.
     * @param[out]  derivatives   The computed derivatives. Must be pre-allocated to the correct size.
     */
    virtual void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                             Eigen::VectorXd &derivatives) const = 0;

    /// Return the refStar error appropriate for this model (e.g. fluxErr or magErr).
    virtual double getRefError(RefStar const &refStar) const = 0;

    /// Return the fittedStar - refStar residual appropriate for this model (e.g. flux - flux or mag - mag).
    virtual double computeRefResidual(FittedStar const &fittedStar, RefStar const &refStar) const = 0;

    /**
     * Return the mapping of ccdImage represented as a PhotoCalib.
     */
    virtual std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const = 0;

    /// Return the number of parameters in the mapping of CcdImage
    unsigned getNpar(CcdImage const &ccdImage) const { return findMapping(ccdImage)->getNpar(); }

    /// Get the mapping associated with ccdImage.
    PhotometryMappingBase const &getMapping(CcdImage const &ccdImage) const {
        return *(findMapping(ccdImage));
    }

    /// Return the total number of parameters in this model.
    virtual int getTotalParameters() const = 0;

    /// Dump the contents of the transforms, for debugging.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /**
     * Return true if this is a "reasonable" model.
     *
     * A valid photometry model is positive within each sensor's bounding box.
     *
     * @param ccdImageList The ccdImages to test the model validity on.
     * @param ndof The number of degrees of freedom in the fit, e.g. from Fitterbase.computeChi2().
     *
     * @return True if the model is valid on all ccdImages.
     */
    bool validate(CcdImageList const &ccdImageList, int ndof) const;

    /**
     * Check that the model is positive on the ccdImage bbox.
     *
     * @param ccdImage The ccdImage to test.
     * @return True if the image is positive on a sampling of points of the ccdImage bbox.
     */
    bool checkPositiveOnBBox(CcdImage const &ccdImage) const;

    friend std::ostream &operator<<(std::ostream &s, PhotometryModel const &model) {
        model.dump(s);
        return s;
    }

    double getErrorPedestal() { return errorPedestal; }

    /// Add a fraction of the instrumental flux to the instrumental flux error, in quadrature.
    double tweakFluxError(jointcal::MeasuredStar const &measuredStar) const {
        if (errorPedestal == 0) {
            return measuredStar.getInstFluxErr();
        } else {
            return std::hypot(measuredStar.getInstFluxErr(), measuredStar.getInstFlux() * errorPedestal);
        }
    }

    /// Add a small magnitude offset to the "instrumental magnitude" error, in quadrature.
    double tweakMagnitudeError(jointcal::MeasuredStar const &measuredStar) const {
        if (errorPedestal == 0) {
            return measuredStar.getInstMagErr();
        } else {
            return std::hypot(measuredStar.getInstMagErr(), errorPedestal);
        }
    }

protected:
    /// Return a pointer to the mapping associated with this ccdImage.
    virtual PhotometryMappingBase *findMapping(CcdImage const &ccdImage) const = 0;

    /// lsst.logging instance, to be created by a subclass so that messages have consistent name.
    LOG_LOGGER _log;

    // Pedestal on flux/magnitude error (percent of flux or delta magnitude)
    double errorPedestal;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MODEL_H
