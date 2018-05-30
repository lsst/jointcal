// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_ASTROMETRY_MODEL_H

#include "memory"

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

    //! Offset the parameters by the provided amounts.
    /*! The shifts are applied according to the indices given in
        AssignIndices. */
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

protected:
    /// Return a pointer to the mapping associated with this ccdImage.
    virtual AstrometryMapping *findMapping(CcdImage const &ccdImage) const = 0;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTROMETRY_MODEL_H
