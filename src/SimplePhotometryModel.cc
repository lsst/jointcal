// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <cmath>

#include "lsst/log/Log.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryTransform.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace lsst {
namespace jointcal {

Eigen::Index SimplePhotometryModel::assignIndices(std::string const &whatToFit, Eigen::Index firstIndex) {
    Eigen::Index ipar = firstIndex;
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
                                              std::vector<Eigen::Index> &indices) const {
    auto mapping = findMapping(ccdImage);
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
}

std::size_t SimplePhotometryModel::getTotalParameters() const {
    std::size_t total = 0;
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

SimpleFluxModel::SimpleFluxModel(CcdImageList const &ccdImageList, double errorPedestal_)
        : SimplePhotometryModel(ccdImageList, LOG_GET("jointcal.SimpleFluxModel"), errorPedestal_) {
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use the single-frame processing calibration from the PhotoCalib as the initial value.
        auto transform = std::make_shared<FluxTransformSpatiallyInvariant>(photoCalib->getCalibrationMean());
        _myMap.emplace(ccdImage->getHashKey(), std::make_unique<PhotometryMapping>(transform));
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
    double tempErr = tweakFluxError(star);
    return mapping->transformError(star, star.getInstFlux(), tempErr);
}

std::shared_ptr<afw::image::PhotoCalib> SimpleFluxModel::toPhotoCalib(CcdImage const &ccdImage) const {
    double calibration = (findMapping(ccdImage)->getParameters()[0]);
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    return std::make_unique<afw::image::PhotoCalib>(calibration, oldPhotoCalib->getCalibrationErr());
}

SimpleMagnitudeModel::SimpleMagnitudeModel(CcdImageList const &ccdImageList, double errorPedestal_)
        : SimplePhotometryModel(ccdImageList, LOG_GET("jointcal.SimpleMagnitudeModel"), errorPedestal_) {
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use the single-frame processing calibration from the PhotoCalib as the default.
        double calib = utils::nanojanskyToABMagnitude(photoCalib->getCalibrationMean());
        auto transform = std::make_shared<MagnitudeTransformSpatiallyInvariant>(calib);
        _myMap.emplace(ccdImage->getHashKey(), std::make_unique<PhotometryMapping>(transform));
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
    double tempErr = tweakMagnitudeError(star);
    return mapping->transformError(star, star.getInstMag(), tempErr);
}

std::shared_ptr<afw::image::PhotoCalib> SimpleMagnitudeModel::toPhotoCalib(CcdImage const &ccdImage) const {
    // NOTE: photocalib is defined as `instFlux * calibration = flux`,
    // so we have to convert the transform from magnitude space.
    double calibration = utils::ABMagnitudeToNanojansky(findMapping(ccdImage)->getParameters()[0]);
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    return std::make_unique<afw::image::PhotoCalib>(calibration, oldPhotoCalib->getCalibrationErr());
}

}  // namespace jointcal
}  // namespace lsst
