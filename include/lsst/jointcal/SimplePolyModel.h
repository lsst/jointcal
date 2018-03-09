// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_POLY_MODEL_H
#define LSST_JOINTCAL_SIMPLE_POLY_MODEL_H

#include "memory"

#include "lsst/jointcal/Eigenstuff.h"

#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include <map>

namespace lsst {
namespace jointcal {

class CcdImage;

/* We deal here with coordinate transforms which are fitted
   and/or necessary to AstrometryFit. The classes SimplePolyModel
and SimplePolyMapping implement a model where  there is one
separate transfrom per CcdImage. One could chose other setups.

*/

//! this is the model used to fit independent CCDs, meaning that there is no instrument model.
/* This modeling of distortions can even accommodate images set mixing instruments */
class SimplePolyModel : public AstrometryModel {
public:
    //! Sky2TP is just a name, it can be anything
    SimplePolyModel(CcdImageList const &ccdImageList, ProjectionHandler const *projectionHandler,
                    bool initFromWcs, unsigned nNotFit = 0, unsigned degree = 3);

    /// No copy or move: there is only ever one instance of a given model (i.e.. per ccd+visit)
    SimplePolyModel(SimplePolyModel const &) = delete;
    SimplePolyModel(SimplePolyModel &&) = delete;
    SimplePolyModel &operator=(SimplePolyModel const &) = delete;
    SimplePolyModel &operator=(SimplePolyModel &&) = delete;

    // The following routines are the interface to AstrometryFit
    //!
    Mapping const *getMapping(CcdImage const &) const;

    //! Positions the various parameter sets into the parameter vector, starting at firstIndex
    unsigned assignIndices(unsigned firstIndex, std::string const &whatToFit);

    // dispaches the offsets after a fit step into the actual locations of parameters
    void offsetParams(Eigen::VectorXd const &delta);

    /*! the mapping of sky coordinates (i.e. the coordinate system
    in which fitted stars are reported) onto the Tangent plane
    (into which the pixel coordinates are transformed) */
    Gtransfo const *getSky2TP(CcdImage const &ccdImage) const { return _sky2TP->getSky2TP(ccdImage); }

    //!
    virtual void freezeErrorTransform();

    //! Access to mappings
    Gtransfo const &getTransfo(CcdImage const &ccdImage) const;

    std::shared_ptr<TanSipPix2RaDec> produceSipWcs(CcdImage const &ccdImage) const;

    ~SimplePolyModel(){};

private:
    typedef std::map<CcdImage const *, std::unique_ptr<SimpleGtransfoMapping>> mapType;
    mapType _myMap;
    ProjectionHandler const *_sky2TP;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_POLY_MODEL_H
