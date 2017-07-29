// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_ASTROMETRY_MODEL_H

#include "memory"

#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/Mapping.h"

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
    //! Mapping associated to a given CcdImage
    virtual const Mapping *getMapping(CcdImage const &) const = 0;

    //! Assign indices to parameters involved in mappings, starting at firstIndex. Returns the highest
    //! assigned index.
    virtual unsigned assignIndices(unsigned firstIndex, std::string const &whatToFit) = 0;

    //! Offset the parameters by the provided amounts.
    /*! The shifts are applied according to the indices given in
        AssignIndices. */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    //! The transformation used to project the positions of FittedStars.
    /*! This defines the coordinate system into which the Mapping of
        this Ccdimage maps the pixel coordinates. */
    virtual const Gtransfo *getSky2TP(CcdImage const &ccdImage) const = 0;

    //! Cook up a SIP WCS.
    virtual std::shared_ptr<TanSipPix2RaDec> produceSipWcs(CcdImage const &ccdImage) const = 0;

    //!
    virtual void freezeErrorScales() = 0;

    virtual ~AstrometryModel(){};
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_ASTROMETRY_MODEL_H
