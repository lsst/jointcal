#include <iostream>
#include <memory>

#include "lsst/log/Log.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.SimplePhotometryModel");
}

namespace lsst {
namespace jointcal {

SimplePhotometryModel::SimplePhotometryModel(CcdImageList const &ccdImageList) {
    _myMap.reserve(ccdImageList.size());
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use the single-frame processing calibration from the PhotoCalib as the default.
        auto transfo =
                std::make_shared<PhotometryTransfoSpatiallyInvariant>(photoCalib->getCalibrationMean());
        _myMap.emplace(ccdImage->getHashKey(), std::make_unique<PhotometryMapping>(transfo));
    }
    LOGLS_INFO(_log, "SimplePhotometryModel got " << _myMap.size() << " ccdImage mappings.");
}

unsigned SimplePhotometryModel::assignIndices(std::string const & /*whatToFit*/, unsigned firstIndex) {
    unsigned ipar = firstIndex;
    for (auto const &i : _myMap) {
        auto mapping = i.second.get();
        mapping->setIndex(ipar);
        ipar += mapping->getNpar();
    }
    return ipar;
}

void SimplePhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

double SimplePhotometryModel::transform(CcdImage const &ccdImage, MeasuredStar const &star,
                                        double instFlux) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(star, instFlux);
}

void SimplePhotometryModel::getMappingIndices(CcdImage const &ccdImage,
                                              std::vector<unsigned> &indices) const {
    auto mapping = findMapping(ccdImage);
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
}

void SimplePhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                        CcdImage const &ccdImage,
                                                        Eigen::VectorXd &derivatives) const {
    auto mapping = findMapping(ccdImage);
    mapping->computeParameterDerivatives(measuredStar, measuredStar.getInstFlux(), derivatives);
}

std::shared_ptr<afw::image::PhotoCalib> SimplePhotometryModel::toPhotoCalib(CcdImage const &ccdImage) const {
    double calibration = (findMapping(ccdImage)->getParameters()[0]);
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    return std::make_unique<afw::image::PhotoCalib>(calibration, oldPhotoCalib->getCalibrationErr());
}

void SimplePhotometryModel::dump(std::ostream &stream) const {
    for (auto &i : _myMap) {
        i.second->dump(stream);
    }
}

PhotometryMappingBase *SimplePhotometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end()) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimplePhotometryModel cannot find CcdImage " + ccdImage.getName());
    }
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
