// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_ASTROMETRY_MODEL_H
#define LSST_JOINTCAL_SIMPLE_ASTROMETRY_MODEL_H

#include "memory"

#include "lsst/jointcal/Eigenstuff.h"

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include <map>

namespace lsst {
namespace jointcal {

class CcdImage;

/* We deal here with coordinate transforms which are fitted
   and/or necessary to AstrometryFit. The classes SimpleAstrometryModel
and SimplePolyMapping implement a model where  there is one
separate transfrom per CcdImage. One could chose other setups.

*/

//! this is the model used to fit independent CCDs, meaning that there is no instrument model.
/* This modeling of distortions can even accommodate images set mixing instruments */
class SimpleAstrometryModel : public AstrometryModel {
public:
    //! Sky2TP is just a name, it can be anything
    SimpleAstrometryModel(CcdImageList const &ccdImageList,
                          const std::shared_ptr<ProjectionHandler const> projectionHandler, bool initFromWCS,
                          unsigned nNotFit = 0, unsigned order = 3);

    /// No copy or move: there is only ever one instance of a given model (i.e.. per ccd+visit)
    SimpleAstrometryModel(SimpleAstrometryModel const &) = delete;
    SimpleAstrometryModel(SimpleAstrometryModel &&) = delete;
    SimpleAstrometryModel &operator=(SimpleAstrometryModel const &) = delete;
    SimpleAstrometryModel &operator=(SimpleAstrometryModel &&) = delete;

    // The following routines are the interface to AstrometryFit
    //!
    const AstrometryMapping *getMapping(CcdImage const &) const;

    //! Positions the various parameter sets into the parameter vector, starting at firstIndex
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex);

    // dispaches the offsets after a fit step into the actual locations of parameters
    void offsetParams(Eigen::VectorXd const &delta);

    /*! the mapping of sky coordinates (i.e. the coordinate system
    in which fitted stars are reported) onto the Tangent plane
    (into which the pixel coordinates are transformed) */
    const std::shared_ptr<Gtransfo const> getSky2TP(CcdImage const &ccdImage) const {
        return _sky2TP->getSky2TP(ccdImage);
    }

    //!
    virtual void freezeErrorTransform();

    /// @copydoc AstrometryModel::getTotalParameters
    int getTotalParameters() const override;

    //! Access to mappings
    Gtransfo const &getTransfo(CcdImage const &ccdImage) const;

    /// @copydoc AstrometryModel::makeSkyWcs
    std::shared_ptr<afw::geom::SkyWcs> makeSkyWcs(CcdImage const &ccdImage) const;

    ~SimpleAstrometryModel(){};

private:
    std::map<CcdImageKey, std::unique_ptr<SimpleGtransfoMapping>> _myMap;
    const std::shared_ptr<ProjectionHandler const> _sky2TP;

    /// @copydoc AstrometryModel::findMapping
    AstrometryMapping *findMapping(CcdImage const &ccdImage) const;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_ASTROMETRY_MODEL_H
