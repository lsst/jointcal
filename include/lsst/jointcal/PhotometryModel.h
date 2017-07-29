// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_PHOTOMETRY_MODEL_H

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
     * Return the "photometric factor" at a given location on a ccdImage.
     *
     * Multiply this by a Calib's flux/magnitude zero-point to get the updated fluxMag0 at that point.
     *
     * @param[in]  ccdImage  The ccdImage to get the photometric factor for.
     * @param[in]  where     Possition on ccdImage in ccd coordinates.
     *
     * @return     The photometric factor at the given location on ccdImage.
     */
    virtual double photomFactor(CcdImage const &ccdImage, Point const &where) const = 0;

    /// Get the mapping associated with ccdImage.
    PhotometryMapping const &getMapping(CcdImage const &ccdImage) const {
        return *(this->findMapping(ccdImage, "getMapping"));
    }

    /**
     * Get how this set of parameters (of length Npar()) map into the "grand" fit.
     *
     * @param[in]     ccdImage  The ccdImage to find the mapping of.
     * @param[out]    indices   The indices of the mapping associated with ccdImage.
     */
    virtual void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) = 0;

    /**
     * Compute the parametric derivatives of this model.
     *
     * @param[in]   measuredStar  The measured star with the position and flux to compute at.
     * @param[in]   ccdImage      The ccdImage containing the measured star, to find the correct mapping.
     * @param[out]  derivatives   The computed derivatives. Must be pre-allocated to the correct size.
     */
    virtual void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                             Eigen::VectorXd &derivatives) = 0;

    // //! number of parameters to be read in indices.size()
    // virtual void setIndicesAndDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
    //                                       std::vector<unsigned> &indices, Eigen::VectorXd &D) = 0;

    virtual ~PhotometryModel(){};

protected:
    /// Return a pointer to the mapping associated with this ccdImage. name is for describing error messages.
    virtual PhotometryMapping *findMapping(CcdImage const &ccdImage, std::string name) const = 0;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MODEL_H
