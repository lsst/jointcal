#include <map>
#include <limits>

#include "lsst/log/Log.h"

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
            // Use (fluxMag0)^-1 from the PhotoCalib as the default.
            auto chipTransfo = std::make_shared<PhotometryTransfoSpatiallyInvariant>(
                    1.0 / photoCalib->getInstFluxMag0());
            _chipMap[chip] = std::unique_ptr<PhotometryMapping>(new PhotometryMapping(chipTransfo));
        }
        // If the visit is not in the map, add it, otherwise continue.
        if (visitPair == _visitMap.end()) {
            auto visitTransfo = std::make_shared<PhotometryTransfoChebyshev>(visitDegree, focalPlaneBBox);
            _visitMap[visit] = std::unique_ptr<PhotometryMapping>(new PhotometryMapping(visitTransfo));
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
                       std::unique_ptr<ChipVisitPhotometryMapping>(
                               new ChipVisitPhotometryMapping(_chipMap[chip], _visitMap[visit])));
    }
    LOGLS_INFO(_log, "Got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                            << " visit mappings; holding chip " << constrainedChip << " fixed.");
    LOGLS_DEBUG(_log, "CcdImage map has " << _myMap.size() << " mappings, with " << _myMap.bucket_count()
                                          << " buckets and a load factor of " << _myMap.load_factor());
}

unsigned ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    // TODO: currently ignoring whatToFit: eventually implement configurability.
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
    auto mapping = this->findMapping(ccdImage, "transform");
    return mapping->transform(measuredStar, instFlux);
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {
    auto mapping = this->findMapping(ccdImage, "getMappingIndices");
    mapping->getMappingIndices(indices);
    // TODO: I think I need a for loop here, from the above value to that +mapping->getNpar()?
}

void ConstrainedPhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                             CcdImage const &ccdImage,
                                                             Eigen::VectorXd &derivatives) {
    auto mapping = this->findMapping(ccdImage, "computeParameterDerivatives");
    mapping->computeParameterDerivatives(measuredStar, measuredStar.getInstFlux(), derivatives);
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

PhotometryMappingBase *ConstrainedPhotometryModel::findMapping(CcdImage const &ccdImage,
                                                               std::string name) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "ConstrainedPhotometryModel::" + name + ", cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
