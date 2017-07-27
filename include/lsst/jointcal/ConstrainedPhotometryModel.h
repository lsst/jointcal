// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H

#include <map>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace lsst {
namespace jointcal {

class ConstrainedPhotometryModel : public PhotometryModel {
public:
    ConstrainedPhotometryModel(CcdImageList const &ccdImageList) {
        for (auto &ccdImage : ccdImageList) {
            _myMap[ccdImage] = std::unique_ptr<PhotometryMapping>(
                    new PhotometryMapping(PhotometryTransfoSpatiallyInvariant()));
        }
    }

    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) { return 0; }

    void offsetParams(Eigen::VectorXd const &delta) {
        for (auto &i : _myMap) {
            i.second->offsetParams(&delta(i.second->getIndex()));
        }
    }

    double photomFactor(CcdImage const &ccdImage, const Point &where) const { return 0; }

    void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {}

    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) {}

private:
    PhotometryMapping *findMapping(CcdImage const &ccdImage, std::string name) const {
        return nullptr;  // waiting on creation of ConstrainedPhotometryModel.cc
        // auto i = _myMap.find(&ccdImage);
        // if (i == _myMap.end())
        //     throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
        //                       "SimplePolyModel::" + name + ", cannot find CcdImage " + ccdImage.getName());
        // return i->second.get();
    }

    typedef std::map<std::shared_ptr<CcdImage>, std::unique_ptr<PhotometryMapping>> mapType;
    mapType _myMap;
    // typedef std::map<VisitIdType, std::unique_ptr<ConstrainedPhotometryMapping>> VisitMapType;
    // VisitMapType _visitMap;
    // typedef std::map<CcdIdType, std::unique_ptr<ConstrainedPhotometryMapping>> ChipMapType;
    // ChipMapType _chipMap;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H
