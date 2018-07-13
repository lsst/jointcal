#include <iostream>
#include <math.h>

#include "lsst/log/Log.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.SimplePhotometryModel");
}

namespace {
double instFluxFromMag(double mag) { return pow(10, mag / -2.5); }
}  // namespace

namespace lsst {
namespace jointcal {

unsigned SimplePhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
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

void SimplePhotometryModel::freezeErrorTransform() {
    for (auto &i : _myMap) {
        i.second->freezeErrorTransform();
    }
}

void SimplePhotometryModel::getMappingIndices(CcdImage const &ccdImage,
                                              std::vector<unsigned> &indices) const {
    auto mapping = findMapping(ccdImage);
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
}

int SimplePhotometryModel::getTotalParameters() const {
    int total = 0;
    for (auto &i : _myMap) {
        total += i.second->getNpar();
    }
    return total;
}

void SimplePhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                        CcdImage const &ccdImage,
                                                        Eigen::VectorXd &derivatives) const {
    auto mapping = findMapping(ccdImage);
    mapping->computeParameterDerivatives(measuredStar, measuredStar.getInstFlux(), derivatives);
}

void SimplePhotometryModel::dump(std::ostream &stream) const {
    for (auto &i : _myMap) {
        stream << i.first << ": ";
        i.second->dump(stream);
        stream << ", ";
    }
}

PhotometryMappingBase *SimplePhotometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimplePhotometryModel cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

SimpleFluxModel::SimpleFluxModel(CcdImageList const &ccdImageList) : SimplePhotometryModel(ccdImageList) {
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use the single-frame processing calibration from the PhotoCalib as the initial value.
        auto transfo = std::make_shared<FluxTransfoSpatiallyInvariant>(photoCalib->getCalibrationMean());
        _myMap.emplace(ccdImage->getHashKey(), std::make_unique<PhotometryMapping>(transfo));
    }
    LOGLS_INFO(_log, "SimpleFluxModel got " << _myMap.size() << " ccdImage mappings.");
}

double SimpleFluxModel::computeResidual(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const {
    return transform(ccdImage, measuredStar) - measuredStar.getFittedStar()->getFlux();
}

double SimpleFluxModel::transform(CcdImage const &ccdImage, MeasuredStar const &star) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(star, star.getInstFlux());
}

double SimpleFluxModel::transformError(CcdImage const &ccdImage, MeasuredStar const &star) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transformError(star, star.getInstFlux(), star.getInstFluxErr());
}

std::shared_ptr<afw::image::PhotoCalib> SimpleFluxModel::toPhotoCalib(CcdImage const &ccdImage) const {
    double calibration = (findMapping(ccdImage)->getParameters()[0]);
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    return std::make_unique<afw::image::PhotoCalib>(calibration, oldPhotoCalib->getCalibrationErr());
}

SimpleMagnitudeModel::SimpleMagnitudeModel(CcdImageList const &ccdImageList)
        : SimplePhotometryModel(ccdImageList) {
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use the single-frame processing calibration from the PhotoCalib as the default.
        double calib = magFromFlux(photoCalib->getCalibrationMean());
        auto transfo = std::make_shared<MagnitudeTransfoSpatiallyInvariant>(calib);
        _myMap.emplace(ccdImage->getHashKey(), std::make_unique<PhotometryMapping>(transfo));
    }
    LOGLS_INFO(_log, "SimpleMagnitudeModel got " << _myMap.size() << " ccdImage mappings.");
}

double SimpleMagnitudeModel::computeResidual(CcdImage const &ccdImage,
                                             MeasuredStar const &measuredStar) const {
    return transform(ccdImage, measuredStar) - measuredStar.getFittedStar()->getMag();
}

double SimpleMagnitudeModel::transform(CcdImage const &ccdImage, MeasuredStar const &star) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(star, star.getInstMag());
}

double SimpleMagnitudeModel::transformError(CcdImage const &ccdImage, MeasuredStar const &star) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transformError(star, star.getInstMag(), star.getInstMagErr());
}

std::shared_ptr<afw::image::PhotoCalib> SimpleMagnitudeModel::toPhotoCalib(CcdImage const &ccdImage) const {
    // NOTE: photocalib is defined as instFlux * calibration = flux
    double calibration = instFluxFromMag(findMapping(ccdImage)->getParameters()[0]);
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    return std::make_unique<afw::image::PhotoCalib>(calibration, oldPhotoCalib->getCalibrationErr());
}

}  // namespace jointcal
}  // namespace lsst
