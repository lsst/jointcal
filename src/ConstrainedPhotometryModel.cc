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

#include <map>
#include <limits>
#include <vector>
#include <string>

#include "lsst/log/Log.h"

#include "astshim.h"
#include "astshim/ChebyMap.h"
#include "lsst/geom.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/math/TransformBoundedField.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ConstrainedPhotometryModel.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace lsst {
namespace jointcal {

Eigen::Index ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit,
                                                       Eigen::Index firstIndex) {
    Eigen::Index index = firstIndex;
    if (whatToFit.find("Model") == std::string::npos) {
        LOGLS_WARN(_log, "assignIndices was called and Model is *not* in whatToFit");
        return index;
    }

    // If we got here, "Model" is definitely in whatToFit.
    _fittingChips = (whatToFit.find("ModelChip") != std::string::npos);
    _fittingVisits = (whatToFit.find("ModelVisit") != std::string::npos);
    // If nothing more than "Model" is specified, it means fit everything.
    if ((!_fittingChips) && (!_fittingVisits)) {
        _fittingChips = _fittingVisits = true;
    }

    if (_fittingChips) {
        for (auto &idMapping : _chipMap) {
            auto mapping = idMapping.second.get();
            // Don't assign indices for fixed parameters.
            if (mapping->isFixed()) continue;
            mapping->setIndex(index);
            index += mapping->getNpar();
        }
    }
    if (_fittingVisits) {
        for (auto &idMapping : _visitMap) {
            auto mapping = idMapping.second.get();
            mapping->setIndex(index);
            index += mapping->getNpar();
        }
    }
    for (auto &idMapping : _chipVisitMap) {
        idMapping.second->setWhatToFit(_fittingChips, _fittingVisits);
    }
    return index;
}

void ConstrainedPhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    if (_fittingChips) {
        for (auto &idMapping : _chipMap) {
            auto mapping = idMapping.second.get();
            // Don't offset indices for fixed parameters.
            if (mapping->isFixed()) continue;
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
    }
    if (_fittingVisits) {
        for (auto &idMapping : _visitMap) {
            auto mapping = idMapping.second.get();
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
    }
}

void ConstrainedPhotometryModel::freezeErrorTransform() {
    for (auto &idMapping : _chipMap) {
        idMapping.second->freezeErrorTransform();
    }
    for (auto &idMapping : _visitMap) {
        idMapping.second->freezeErrorTransform();
    }
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage, IndexVector &indices) const {
    auto mapping = findMapping(ccdImage);
    mapping->getMappingIndices(indices);
}

std::size_t ConstrainedPhotometryModel::getTotalParameters() const {
    std::size_t total = 0;
    for (auto &idMapping : _chipMap) {
        total += idMapping.second->getNpar();
    }
    for (auto &idMapping : _visitMap) {
        total += idMapping.second->getNpar();
    }
    return total;
}

void ConstrainedPhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                             CcdImage const &ccdImage,
                                                             Eigen::VectorXd &derivatives) const {
    auto mapping = findMapping(ccdImage);
    mapping->computeParameterDerivatives(measuredStar, measuredStar.getInstFlux(), derivatives);
}

namespace {
// Convert photoTransform's way of storing Chebyshev coefficients into the format wanted by ChebyMap.
ndarray::Array<double, 2, 2> toChebyMapCoeffs(const std::shared_ptr<PhotometryTransformChebyshev>& transform) {
    auto coeffs = transform->getCoefficients();
    // 4 x nPar: ChebyMap wants rows that look like (a_ij, 1, i, j) for out += a_ij*T_i(x)*T_j(y)
    ndarray::Array<double, 2, 2> chebyCoeffs =
            allocate(ndarray::makeVector(transform->getNpar(), std::size_t(4)));
    Eigen::VectorXd::Index k = 0;
    auto order = transform->getOrder();
    for (ndarray::Size j = 0; j <= order; ++j) {
        ndarray::Size const iMax = order - j;  // to save re-computing `i+j <= order` every inner step.
        for (ndarray::Size i = 0; i <= iMax; ++i, ++k) {
            chebyCoeffs[k][0] = coeffs[j][i];
            chebyCoeffs[k][1] = 1;
            chebyCoeffs[k][2] = i;
            chebyCoeffs[k][3] = j;
        }
    }
    return chebyCoeffs;
}
}  // namespace

void ConstrainedPhotometryModel::print(std::ostream &out) const {
    for (auto &idMapping : _chipMap) {
        out << "Sensor: " << idMapping.first << std::endl;
        idMapping.second->print(out);
        out << std::endl;
    }
    out << std::endl;
    for (auto &idMapping : _visitMap) {
        out << "Visit: " << idMapping.first << std::endl;
        idMapping.second->print(out);
        out << std::endl;
    }
}

PhotometryMappingBase *ConstrainedPhotometryModel::findMapping(CcdImage const &ccdImage) const {
    auto idMapping = _chipVisitMap.find(ccdImage.getHashKey());
    if (idMapping == _chipVisitMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "ConstrainedPhotometryModel cannot find CcdImage " + ccdImage.getName());
    return idMapping->second.get();
}

template <class ChipTransform, class VisitTransform, class ChipVisitMapping>
void ConstrainedPhotometryModel::initialize(CcdImageList const &ccdImageList,
                                            geom::Box2D const &focalPlaneBBox, int visitOrder) {
    // keep track of which chip we want to constrain (the one closest to the middle of the focal plane)
    double minRadius2 = std::numeric_limits<double>::infinity();
    CcdIdType constrainedChip = -1;

    // First initialize all visit and ccd transforms, before we make the ccdImage mappings.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        auto visitPair = _visitMap.find(visit);
        auto chipPair = _chipMap.find(chip);

        // If the chip is not in the map, add it, otherwise continue.
        if (chipPair == _chipMap.end()) {
            auto center = ccdImage->getDetector()->getCenter(afw::cameraGeom::FOCAL_PLANE);
            double radius2 = std::pow(center.getX(), 2) + std::pow(center.getY(), 2);
            if (radius2 < minRadius2) {
                minRadius2 = radius2;
                constrainedChip = chip;
            }
            auto photoCalib = ccdImage->getPhotoCalib();
            // Use the single-frame processing calibration from the PhotoCalib as the default.
            auto chipTransform = std::make_unique<ChipTransform>(initialChipCalibration(photoCalib));
            _chipMap[chip] = std::make_shared<PhotometryMapping>(std::move(chipTransform));
        }
        // If the visit is not in the map, add it, otherwise continue.
        if (visitPair == _visitMap.end()) {
            auto visitTransform = std::make_unique<VisitTransform>(visitOrder, focalPlaneBBox);
            _visitMap[visit] = std::make_shared<PhotometryMapping>(std::move(visitTransform));
        }
    }

    // Fix one chip mapping, to remove the degeneracy from the system.
    _chipMap.at(constrainedChip)->setFixed(true);

    // Now create the ccdImage mappings, which are combinations of the chip/visit mappings above.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        _chipVisitMap.emplace(ccdImage->getHashKey(),
                              std::make_unique<ChipVisitMapping>(_chipMap[chip], _visitMap[visit]));
    }
    LOGLS_INFO(_log, "Got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                            << " visit mappings; holding chip " << constrainedChip << " fixed ("
                            << getTotalParameters() << " total parameters).");
    LOGLS_DEBUG(_log, "CcdImage map has " << _chipVisitMap.size() << " mappings, with "
                                          << _chipVisitMap.bucket_count() << " buckets and a load factor of "
                                          << _chipVisitMap.load_factor());
}

ConstrainedPhotometryModel::PrepPhotoCalib ConstrainedPhotometryModel::prepPhotoCalib(
        CcdImage const &ccdImage) const {
    auto detector = ccdImage.getDetector();
    auto ccdBBox = detector->getBBox();
    auto *mapping = dynamic_cast<ChipVisitPhotometryMapping *>(findMapping(ccdImage));

    // There should be no way in which we can get to this point and not have a ChipVisitMapping,
    // so blow up if we don't.
    assert(mapping != nullptr);
    // We know it's a Chebyshev transform because we created it as such, so blow up if it's not.
    auto visitPhotometryTransform = std::dynamic_pointer_cast<PhotometryTransformChebyshev>(
            mapping->getVisitMapping()->getTransform());
    assert(visitPhotometryTransform != nullptr);
    auto focalBBox = visitPhotometryTransform->getBBox();

    // Unravel our chebyshev coefficients to build an astshim::ChebyMap.
    auto coeff_f = toChebyMapCoeffs(std::dynamic_pointer_cast<PhotometryTransformChebyshev>(
            mapping->getVisitMapping()->getTransform()));
    // Bounds are the bbox
    std::vector<double> lowerBound = {focalBBox.getMinX(), focalBBox.getMinY()};
    std::vector<double> upperBound = {focalBBox.getMaxX(), focalBBox.getMaxY()};
    afw::geom::TransformPoint2ToGeneric visitTransform(ast::ChebyMap(coeff_f, 1, lowerBound, upperBound));

    double chipConstant = mapping->getChipMapping()->getParameters()[0];

    // Compute a box that covers the area of the ccd in focal plane coordinates.
    // This is the box over which we want to compute the mean of the visit transform.
    auto pixToFocal = detector->getTransform(afw::cameraGeom::PIXELS, afw::cameraGeom::FOCAL_PLANE);
    geom::Box2D ccdBBoxInFocal;
    for (auto const &point : pixToFocal->applyForward(geom::Box2D(ccdBBox).getCorners())) {
        ccdBBoxInFocal.include(point);
    }
    double visitMean = visitPhotometryTransform->mean(ccdBBoxInFocal);

    return {chipConstant, visitTransform, pixToFocal, visitMean};
}

// ConstrainedFluxModel methods

double ConstrainedFluxModel::computeResidual(CcdImage const &ccdImage,
                                             MeasuredStar const &measuredStar) const {
    return transform(ccdImage, measuredStar) - measuredStar.getFittedStar()->getFlux();
}

double ConstrainedFluxModel::transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(measuredStar, measuredStar.getInstFlux());
}

double ConstrainedFluxModel::transformError(CcdImage const &ccdImage,
                                            MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    double tempErr = tweakFluxError(measuredStar);
    return mapping->transformError(measuredStar, measuredStar.getInstFlux(), tempErr);
}

std::shared_ptr<afw::image::PhotoCalib> ConstrainedFluxModel::toPhotoCalib(CcdImage const &ccdImage) const {
    auto ccdBBox = ccdImage.getDetector()->getBBox();
    auto prep = prepPhotoCalib(ccdImage);

    // The chip part is easy: zoom map with the single value as the "zoom" factor
    afw::geom::Transform<afw::geom::GenericEndpoint, afw::geom::GenericEndpoint> zoomTransform(
            ast::ZoomMap(1, prep.chipConstant));

    // Now stitch them all together.
    auto transform = prep.pixToFocal->then(prep.visitTransform)->then(zoomTransform);

    // NOTE: TransformBoundedField does not implement mean(), so we have to compute it here.
    double mean = prep.chipConstant * prep.visitMean;

    auto boundedField = std::make_shared<afw::math::TransformBoundedField>(ccdBBox, *transform);
    return std::make_shared<afw::image::PhotoCalib>(mean, ccdImage.getPhotoCalib()->getCalibrationErr(),
                                                    boundedField, false);
}

void ConstrainedFluxModel::print(std::ostream &out) const {
    out << "ConstrainedFluxModel:" << std::endl;
    ConstrainedPhotometryModel::print(out);
}

// ConstrainedMagnitudeModel methods

double ConstrainedMagnitudeModel::computeResidual(CcdImage const &ccdImage,
                                                  MeasuredStar const &measuredStar) const {
    return transform(ccdImage, measuredStar) - measuredStar.getFittedStar()->getMag();
}

double ConstrainedMagnitudeModel::transform(CcdImage const &ccdImage,
                                            MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(measuredStar, measuredStar.getInstMag());
}

double ConstrainedMagnitudeModel::transformError(CcdImage const &ccdImage,
                                                 MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    double tempErr = tweakFluxError(measuredStar);
    return mapping->transformError(measuredStar, measuredStar.getInstFlux(), tempErr);
}

std::shared_ptr<afw::image::PhotoCalib> ConstrainedMagnitudeModel::toPhotoCalib(
        CcdImage const &ccdImage) const {
    auto ccdBBox = ccdImage.getDetector()->getBBox();
    auto prep = prepPhotoCalib(ccdImage);

    using namespace std::string_literals;  // for operator""s to convert string literal->std::string
    afw::geom::Transform<afw::geom::GenericEndpoint, afw::geom::GenericEndpoint> logTransform(
            ast::MathMap(1, 1, {"y=pow(10.0,x/-2.5)"s}, {"x=-2.5*log10(y)"s}));

    // The chip part is easy: zoom map with the value (converted to a flux) as the "zoom" factor.
    double chipCalibration = utils::ABMagnitudeToNanojansky(prep.chipConstant);
    afw::geom::Transform<afw::geom::GenericEndpoint, afw::geom::GenericEndpoint> zoomTransform(
            ast::ZoomMap(1, chipCalibration));

    // Now stitch them all together.
    auto transform = prep.pixToFocal->then(prep.visitTransform)->then(logTransform)->then(zoomTransform);

    // NOTE: TransformBoundedField does not implement mean(), so we have to compute it here.
    double mean = chipCalibration * std::pow(10, prep.visitMean / -2.5);

    auto boundedField = std::make_shared<afw::math::TransformBoundedField>(ccdBBox, *transform);
    return std::make_shared<afw::image::PhotoCalib>(mean, ccdImage.getPhotoCalib()->getCalibrationErr(),
                                                    boundedField, false);
}

void ConstrainedMagnitudeModel::print(std::ostream &out) const {
    out << "ConstrainedMagnitudeModel (" << _chipVisitMap.size() << " composite mappings; " << _chipMap.size()
        << " sensor mappings, " << _visitMap.size() << " visit mappings):" << std::endl;
    ConstrainedPhotometryModel::print(out);
}

// explicit instantiation of templated function, so pybind11 can
template void ConstrainedPhotometryModel::initialize<FluxTransformSpatiallyInvariant, FluxTransformChebyshev,
                                                     ChipVisitFluxMapping>(CcdImageList const &,
                                                                           geom::Box2D const &, int);
template void ConstrainedPhotometryModel::initialize<MagnitudeTransformSpatiallyInvariant,
                                                     MagnitudeTransformChebyshev, ChipVisitMagnitudeMapping>(
        CcdImageList const &, geom::Box2D const &, int);

}  // namespace jointcal
}  // namespace lsst
