#include <map>
#include <limits>

#include "lsst/log/Log.h"

#include "astshim.h"
#include "astshim/ChebyMap.h"
#include "lsst/afw/geom/Transform.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/math/TransformBoundedField.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ConstrainedPhotometryModel.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.ConstrainedPhotometryModel");
}

namespace lsst {
namespace jointcal {

ConstrainedPhotometryModel::ConstrainedPhotometryModel(CcdImageList const &ccdImageList,
                                                       afw::geom::Box2D const &focalPlaneBBox, int visitOrder)
        : _fittingChips(false), _fittingVisits(false) {
    // keep track of which chip we want to constrain (the one closest to the middle of the focal plane)
    double minRadius2 = std::numeric_limits<double>::infinity();
    CcdIdType constrainedChip = -1;

    // First initialize all visit and ccd transfos, before we make the ccdImage mappings.
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
            auto chipTransfo =
                    std::make_shared<FluxTransfoSpatiallyInvariant>(photoCalib->getCalibrationMean());
            _chipMap[chip] = std::make_unique<PhotometryMapping>(std::move(chipTransfo));
        }
        // If the visit is not in the map, add it, otherwise continue.
        if (visitPair == _visitMap.end()) {
            auto visitTransfo = std::make_shared<PhotometryTransfoChebyshev>(visitOrder, focalPlaneBBox);
            _visitMap[visit] = std::make_unique<PhotometryMapping>(std::move(visitTransfo));
        }
    }

    // Fix one chip mapping, to remove the degeneracy from the system.
    _chipMap.at(constrainedChip)->setFixed(true);

    // Now create the ccdImage mappings, which are combinations of the chip/visit mappings above.
    _chipVisitMap.reserve(
            ccdImageList.size());  // we know how chipVisitbig it will be, so pre-allocate space.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        _chipVisitMap.emplace(ccdImage->getHashKey(),
                              std::make_unique<ChipVisitPhotometryMapping>(_chipMap[chip], _visitMap[visit]));
    }
    LOGLS_INFO(_log, "Got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                            << " visit mappings; holding chip " << constrainedChip << " fixed ("
                            << getTotalParameters() << " total parameters).");
    LOGLS_DEBUG(_log, "CcdImage map has " << _chipVisitMap.size() << " mappings, with "
                                          << _chipVisitMap.bucket_count() << " buckets and a load factor of "
                                          << _chipVisitMap.load_factor());
}

unsigned ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    unsigned index = firstIndex;
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

double ConstrainedPhotometryModel::computeResidual(CcdImage const &ccdImage,
                                                   MeasuredStar const &measuredStar) const {
    return transform(ccdImage, measuredStar) - measuredStar.getFittedStar()->getFlux();
}

double ConstrainedPhotometryModel::transform(CcdImage const &ccdImage,
                                             MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(measuredStar, measuredStar.getInstFlux());
}

double ConstrainedPhotometryModel::transformError(CcdImage const &ccdImage,
                                                  MeasuredStar const &measuredStar) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transformError(measuredStar, measuredStar.getInstFluxErr());
}

void ConstrainedPhotometryModel::freezeErrorTransform() {
    for (auto &idMapping : _chipMap) {
        idMapping.second.get()->freezeErrorTransform();
    }
    for (auto &idMapping : _visitMap) {
        idMapping.second.get()->freezeErrorTransform();
    }
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage,
                                                   std::vector<unsigned> &indices) const {
    auto mapping = findMapping(ccdImage);
    mapping->getMappingIndices(indices);
}

int ConstrainedPhotometryModel::getTotalParameters() const {
    int total = 0;
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
// Convert photoTransfo's way of storing Chebyshev coefficients into the format wanted by ChebyMap.
ndarray::Array<double, 2, 2> toChebyMapCoeffs(std::shared_ptr<PhotometryTransfoChebyshev> transfo) {
    auto coeffs = transfo->getCoefficients();
    // 4 x nPar: ChebyMap wants rows that look like (a_ij, 1, i, j) for out += a_ij*T_i(x)*T_j(y)
    ndarray::Array<double, 2, 2> chebyCoeffs = allocate(ndarray::makeVector(transfo->getNpar(), 4));
    Eigen::VectorXd::Index k = 0;
    auto order = transfo->getOrder();
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

std::shared_ptr<afw::image::PhotoCalib> ConstrainedPhotometryModel::toPhotoCalib(
        CcdImage const &ccdImage) const {
    auto oldPhotoCalib = ccdImage.getPhotoCalib();
    auto detector = ccdImage.getDetector();
    auto ccdBBox = detector->getBBox();
    ChipVisitPhotometryMapping *mapping = dynamic_cast<ChipVisitPhotometryMapping *>(findMapping(ccdImage));
    // There should be no way in which we can get to this point and not have a ChipVisitMapping,
    // so blow up if we don't.
    assert(mapping != nullptr);
    auto pixToFocal = detector->getTransform(afw::cameraGeom::PIXELS, afw::cameraGeom::FOCAL_PLANE);
    // We know it's a Chebyshev transfo because we created it as such, so blow up if it's not.
    auto visitTransfo =
            std::dynamic_pointer_cast<PhotometryTransfoChebyshev>(mapping->getVisitMapping()->getTransfo());
    assert(visitTransfo != nullptr);
    auto focalBBox = visitTransfo->getBBox();

    // Unravel our chebyshev coefficients to build an astshim::ChebyMap.
    auto coeff_f = toChebyMapCoeffs(
            std::dynamic_pointer_cast<PhotometryTransfoChebyshev>(mapping->getVisitMapping()->getTransfo()));
    // Bounds are the bbox
    std::vector<double> lowerBound = {focalBBox.getMinX(), focalBBox.getMinY()};
    std::vector<double> upperBound = {focalBBox.getMaxX(), focalBBox.getMaxY()};

    afw::geom::TransformPoint2ToGeneric chebyTransform(ast::ChebyMap(coeff_f, 1, lowerBound, upperBound));

    // The chip part is easy: zoom map with the single value as the "zoom" factor.
    afw::geom::Transform<afw::geom::GenericEndpoint, afw::geom::GenericEndpoint> zoomTransform(
            ast::ZoomMap(1, mapping->getChipMapping()->getParameters()[0]));

    // Now stitch them all together.
    auto transform = pixToFocal->then(chebyTransform)->then(zoomTransform);
    // NOTE: TransformBoundedField does not yet implement mean(), so we have to compute it here.
    double mean = mapping->getChipMapping()->getParameters()[0] * visitTransfo->mean();
    auto boundedField = std::make_shared<afw::math::TransformBoundedField>(ccdBBox, *transform);
    return std::make_shared<afw::image::PhotoCalib>(mean, oldPhotoCalib->getCalibrationErr(), boundedField,
                                                    false);
}

void ConstrainedPhotometryModel::dump(std::ostream &stream) const {
    for (auto &idMapping : _chipMap) {
        idMapping.second->dump(stream);
        stream << std::endl;
    }
    stream << std::endl;
    for (auto &idMapping : _visitMap) {
        idMapping.second->dump(stream);
        stream << std::endl;
    }
}

PhotometryMappingBase *ConstrainedPhotometryModel::findMapping(CcdImage const &ccdImage) const {
    auto idMapping = _chipVisitMap.find(ccdImage.getHashKey());
    if (idMapping == _chipVisitMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "ConstrainedPhotometryModel cannot find CcdImage " + ccdImage.getName());
    return idMapping->second.get();
}

}  // namespace jointcal
}  // namespace lsst
