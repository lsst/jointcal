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
    ConstrainedPhotometryModel(const CcdImageList &ccdImageList);

    unsigned assignIndices(const std::string &whatToFit, unsigned firstIndex);

    void offsetParams(const Eigen::VectorXd &delta);

    double photomFactor(const CcdImage& ccdImage, const Point &where) const;

    void getIndicesAndDerivatives(const MeasuredStar &measuredStar,
                                  const CcdImage &ccdImage,
                                  std::vector<unsigned> &indices,
                                  Eigen::VectorXd &D);

private:
    typedef std::map<std::shared_ptr<CcdImage>, std::unique_ptr<PhotometryMapping>> mapType;
    // typedef std::map<VisitIdType, std::unique_ptr<ConstrainedPhotometryMapping>> VisitMapType;
    // VisitMapType _visitMap;
    // typedef std::map<CcdIdType, std::unique_ptr<ConstrainedPhotometryMapping>> ChipMapType;
    // ChipMapType _chipMap;
};

}} // namespaces

#endif // LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H
