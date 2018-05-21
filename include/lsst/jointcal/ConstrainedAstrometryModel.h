// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H

#include "memory"  // for std::*_ptr

#include "lsst/jointcal/Eigenstuff.h"

class CcdImage;

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"
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
    AstrometryMapping const *getMapping(CcdImage const &) const;

    /**
     * Positions the various parameter sets into the parameter vector, starting at
     * firstIndex.
     */
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex);

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

    /// @copydoc AstrometryModel::getTotalParameters
    int getTotalParameters() const override;

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
    const std::shared_ptr<Gtransfo const> getSky2TP(CcdImage const &ccdImage) const {
        return _sky2TP->getSky2TP(ccdImage);
    }

    /// @copydoc AstrometryModel::makeSkyWcs
    std::shared_ptr<afw::geom::SkyWcs> makeSkyWcs(CcdImage const &ccdImage) const;

private:
    std::unordered_map<CcdImageKey, std::unique_ptr<TwoTransfoMapping>> _mappings;
    std::map<CcdIdType, std::shared_ptr<SimpleGtransfoMapping>> _chipMap;
    std::map<VisitIdType, std::shared_ptr<SimpleGtransfoMapping>> _visitMap;
    const std::shared_ptr<ProjectionHandler const> _sky2TP;
    bool _fittingChips, _fittingVisits;

    /// @copydoc AstrometryModel::findMapping
    AstrometryMapping *findMapping(CcdImage const &ccdImage) const;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CONSTRAINED_ASTROMETRY_MODEL_H
