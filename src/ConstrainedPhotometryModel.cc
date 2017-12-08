#include <map>
#include <limits>
#include <memory>

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
                                                       afw::geom::Box2D const &focalPlaneBBox,
                                                       int visitDegree) {
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
            double radius2 = std::pow(center.getPoint().getX(), 2) + std::pow(center.getPoint().getY(), 2);
            if (radius2 < minRadius2) {
                minRadius2 = radius2;
                constrainedChip = chip;
            }
            auto photoCalib = ccdImage->getPhotoCalib();
            // Use the single-frame processing calibration from the PhotoCalib as the default.
            auto chipTransfo =
                    std::make_shared<PhotometryTransfoSpatiallyInvariant>(photoCalib->getCalibrationMean());
            _chipMap[chip] =
                    std::make_unique<PhotometryMapping>(std::move(chipTransfo));
        }
        // If the visit is not in the map, add it, otherwise continue.
        if (visitPair == _visitMap.end()) {
            auto visitTransfo = std::make_shared<PhotometryTransfoChebyshev>(visitDegree, focalPlaneBBox);
            _visitMap[visit] =
                    std::make_unique<PhotometryMapping>(std::move(visitTransfo));
        }
    }

    // Fix one chip mapping, to remove the degeneracy from the system.
    _chipMap.at(constrainedChip)->setFixed(true);

    // Now create the ccdImage mappings, which are combinations of the chip/visit mappings above.
    _myMap.reserve(ccdImageList.size());  // we know how big it will be, so pre-allocate space.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        _myMap.emplace(ccdImage->getHashKey(),
                       std::make_unique<ChipVisitPhotometryMapping>(
                               _chipMap[chip], _visitMap[visit]));
    }
    LOGLS_INFO(_log, "Got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                            << " visit mappings; holding chip " << constrainedChip << " fixed.");
    LOGLS_DEBUG(_log, "CcdImage map has " << _myMap.size() << " mappings, with " << _myMap.bucket_count()
                                          << " buckets and a load factor of " << _myMap.load_factor());
}

unsigned ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    // TODO DM-8046: currently ignoring whatToFit: eventually implement configurability.
    unsigned index = firstIndex;
    for (auto &i : _chipMap) {
        auto mapping = i.second.get();
        // Don't assign indices for fixed parameters.
        if (mapping->isFixed()) continue;
        mapping->setIndex(index);
        index += mapping->getNpar();
    }
    for (auto &i : _visitMap) {
        auto mapping = i.second.get();
        mapping->setIndex(index);
        index += mapping->getNpar();
    }
    return index;
}

void ConstrainedPhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _chipMap) {
        auto mapping = i.second.get();
        // Don't offset indices for fixed parameters.
        if (mapping->isFixed()) continue;
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
    for (auto &i : _visitMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

double ConstrainedPhotometryModel::transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar,
                                             double instFlux) const {
    auto mapping = findMapping(ccdImage);
    return mapping->transform(measuredStar, instFlux);
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage,
                                                   std::vector<unsigned> &indices) const {
    auto mapping = findMapping(ccdImage);
    mapping->getMappingIndices(indices);
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
    auto degree = transfo->getDegree();
    for (ndarray::Size j = 0; j <= degree; ++j) {
        ndarray::Size const iMax = degree - j;  // to save re-computing `i+j <= degree` every inner step.
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
    auto *mapping = dynamic_cast<ChipVisitPhotometryMapping *>(findMapping(ccdImage));
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
    ast::ChebyMap chebyMap(coeff_f, 1, lowerBound, upperBound);

    // The chip part is easy: zoom map with the single value as the "zoom" factor.
    ast::ZoomMap zoomMap(1, mapping->getChipMapping()->getParameters()[0]);

    // Now stitch them all together.
    auto transform =
            afw::geom::TransformPoint2ToGeneric(pixToFocal->getFrameSet()->then(chebyMap).then(zoomMap));
    // NOTE: TransformBoundedField does not yet implement mean(), so we have to compute it here.
    double mean = mapping->getChipMapping()->getParameters()[0] * visitTransfo->mean();
    auto boundedField = std::make_shared<afw::math::TransformBoundedField>(ccdBBox, transform);
    return std::make_shared<afw::image::PhotoCalib>(mean, oldPhotoCalib->getCalibrationErr(), boundedField,
                                                    false);
}

void ConstrainedPhotometryModel::dump(std::ostream &stream) const {
    for (auto &i : _chipMap) {
        i.second->dump(stream);
        stream << std::endl;
    }
    stream << std::endl;
    for (auto &i : _visitMap) {
        i.second->dump(stream);
        stream << std::endl;
    }
}

PhotometryMappingBase *ConstrainedPhotometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "ConstrainedPhotometryModel cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
