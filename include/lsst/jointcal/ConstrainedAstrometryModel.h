// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CONSTRAINED_POLY_MODEL_H
#define LSST_JOINTCAL_CONSTRAINED_POLY_MODEL_H

#include "memory"  // for std::*_ptr

#include "lsst/jointcal/Eigenstuff.h"

class CcdImage;

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/TwoTransfoMapping.h"
#include "lsst/jointcal/CcdImage.h"

#include <map>

namespace lsst {
namespace jointcal {

/**
 * This is the model used to fit mappings as the combination of a
 * transformation depending on the chip number (instrument model) and a
 * transformation per visit (anamorphism). The two-transformation Mapping
 * required for this model is TwoTransfoMapping. This modeling of distortions
 * is meant for a set of images from a single mosaic imager.
 */
class ConstrainedAstrometryModel : public AstrometryModel {
public:
    ConstrainedAstrometryModel(CcdImageList const &ccdImageList, ProjectionHandler const *projectionHandler,
                               bool initFromWCS, int chipDegree = 3, int visitDegree = 2);

    /// No copy or move: there is only ever one instance of a given model (i.e. per ccd+visit)
    ConstrainedAstrometryModel(ConstrainedAstrometryModel const &) = delete;
    ConstrainedAstrometryModel(ConstrainedAstrometryModel &&) = delete;
    ConstrainedAstrometryModel &operator=(ConstrainedAstrometryModel const &) = delete;
    ConstrainedAstrometryModel &operator=(ConstrainedAstrometryModel &&) = delete;

    // The following routines are the interface to AstrometryFit
    //!
    Mapping const *getMapping(CcdImage const &) const;

    /**
     * Positions the various parameter sets into the parameter vector, starting at
     * firstIndex.
     */
    unsigned assignIndices(unsigned firstIndex, std::string const &whatToFit);

    /**
     * Dispaches the offsets after a fit step into the actual locations of
     * parameters.
     */
    void offsetParams(Eigen::VectorXd const &Delta);

    /**
     * From there on, measurement errors are propagated using the current
     * transfos (and no longer evolve).
     */
    void freezeErrorTransform();

    //! Access to mappings
    Gtransfo const &getChipTransfo(CcdIdType const chip) const;

    //! Access to mappings
    Gtransfo const &getVisitTransfo(VisitIdType const &visit) const;

    //! Access to array of visits involved in the solution.
    std::vector<VisitIdType> getVisits() const;

    /**
     * The mapping of sky coordinates (i.e. the coordinate system in which fitted
     * stars are reported) onto the Tangent plane (into which the pixel coordinates
     * are transformed).
     */
    const Gtransfo *getSky2TP(CcdImage const &ccdImage) const { return _sky2TP->getSky2TP(ccdImage); }

    std::shared_ptr<TanSipPix2RaDec> produceSipWcs(CcdImage const &ccdImage) const;

private:
    typedef std::map<const CcdImage *, std::unique_ptr<TwoTransfoMapping>> mappingMapType;
    mappingMapType _mappings;
    typedef std::map<CcdIdType, std::unique_ptr<SimpleGtransfoMapping>> chipMapType;
    chipMapType _chipMap;
    typedef std::map<VisitIdType, std::unique_ptr<SimpleGtransfoMapping>> visitMapType;
    visitMapType _visitMap;
    const ProjectionHandler *_sky2TP;
    bool _fittingChips, _fittingVisits;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CONSTRAINED_POLY_MODEL_H
