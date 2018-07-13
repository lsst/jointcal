// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_PHOTOMETRY_MODEL_H

#include "lsst/afw/image/PhotoCalib.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include <string>
#include <vector>

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;
class MeasuredStar;

//! Interface class for PhotometryFit
class PhotometryModel {
public:
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
     * Offset the parameters by the provided amounts.
     *
     * The shifts are applied according to the indices given in assignIndices.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

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

    /// Dump the contents of the transfos, for debugging.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    friend std::ostream &operator<<(std::ostream &s, PhotometryModel const &model) {
        model.dump(s);
        return s;
    }

protected:
    /// Return a pointer to the mapping associated with this ccdImage.
    virtual PhotometryMappingBase *findMapping(CcdImage const &ccdImage) const = 0;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MODEL_H
