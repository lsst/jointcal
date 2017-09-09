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
     * Assign indices to parameters involved in mappings, starting at firstIndex.
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

    /**
     * Return the on-sky transformed flux for measuredStar on ccdImage.
     *
     * @param[in]  ccdImage     The ccdImage where measuredStar resides.
     * @param      measuredStar The measured star position to compute the transform at.
     * @param[in]  instFlux     The instrument flux to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
     */
    virtual double transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar,
                             double instFlux) const = 0;

    /**
     * Get how this set of parameters (of length Npar()) map into the "grand" fit.
     *
     * @param[out]    indices   The indices of the mapping associated with ccdImage.
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
    // {
    //     return this->findMapping(ccdImage, "getMapping")->toPhotoCalib();
    // }

    /// Return the number of parameters in the mapping of CcdImage
    unsigned getNpar(CcdImage const &ccdImage) const {
        return this->findMapping(ccdImage, "getNpar")->getNpar();
    }

    /// Get the mapping associated with ccdImage.
    PhotometryMappingBase const &getMapping(CcdImage const &ccdImage) const {
        return *(this->findMapping(ccdImage, "getMapping"));
    }

    /// Dump the contents of the transfos, for debugging.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /// For python print().
    std::string __str__() const {
        std::stringstream s;
        dump(s);
        return s.str();
    }

protected:
    /// Return a pointer to the mapping associated with this ccdImage. name is for describing error messages.
    virtual PhotometryMappingBase *findMapping(CcdImage const &ccdImage, std::string name) const = 0;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MODEL_H
