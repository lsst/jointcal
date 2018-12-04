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

#include "astshim.h"
#include "lsst/afw/geom.h"
#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/ConstrainedAstrometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/StarMatch.h"

#include "lsst/pex/exceptions.h"
namespace pexExcept = lsst::pex::exceptions;

#include <memory>
#include <string>
#include <iostream>

namespace {
// Append the keys of this map into a comma-separated string.
template <typename KeyType, typename ValueType>
void outputMapKeys(std::map<KeyType, ValueType> const &map, std::ostream &os) {
    bool first = true;
    os << "[";
    for (auto const &i : map) {
        if (first)
            first = false;
        else
            os << ", ";
        os << i.first;
    }
    os << "]";
}
}  // namespace

namespace lsst {
namespace jointcal {

ConstrainedAstrometryModel::ConstrainedAstrometryModel(
        CcdImageList const &ccdImageList, std::shared_ptr<ProjectionHandler const> projectionHandler,
        int chipOrder, int visitOrder)
        : AstrometryModel(LOG_GET("jointcal.ConstrainedAstrometryModel")),
          _skyToTangentPlane(projectionHandler) {
    // keep track of which chip we want to hold fixed (the one closest to the middle of the focal plane)
    double minRadius2 = std::numeric_limits<double>::infinity();
    CcdIdType constrainedChip = -1;

    // first loop to initialize all visit and chip transforms.
    for (auto &ccdImage : ccdImageList) {
        const CcdImage &im = *ccdImage;
        auto visit = im.getVisit();
        auto chip = im.getCcdId();
        auto visitp = _visitMap.find(visit);
        if (visitp == _visitMap.end()) {
            _visitMap[visit] = std::make_shared<SimplePolyMapping>(AstrometryTransformLinear(),
                                                                   AstrometryTransformPolynomial(visitOrder));
        }
        auto chipp = _chipMap.find(chip);
        if (chipp == _chipMap.end()) {
            auto center = ccdImage->getDetector()->getCenter(afw::cameraGeom::FOCAL_PLANE);
            double radius2 = std::pow(center.getX(), 2) + std::pow(center.getY(), 2);
            if (radius2 < minRadius2) {
                minRadius2 = radius2;
                constrainedChip = chip;
            }

            auto pixelsToFocal =
                    im.getDetector()->getTransform(afw::cameraGeom::PIXELS, afw::cameraGeom::FOCAL_PLANE);
            Frame const &frame = im.getImageFrame();
            // construct the chip transform by approximating the pixel->Focal afw::geom::Transform.
            AstrometryTransformPolynomial pol =
                    AstrometryTransformPolynomial(pixelsToFocal, frame, chipOrder);
            AstrometryTransformLinear shiftAndNormalize = normalizeCoordinatesTransform(frame);
            _chipMap[chip] = std::make_shared<SimplePolyMapping>(shiftAndNormalize,
                                                                 pol * shiftAndNormalize.inverted());
        }
    }

    // Hold the "central" chip map fixed and don't fit it, to remove a degeneracy.
    _chipMap.at(constrainedChip)->setToBeFit(false);

    // now, second loop to set the mappings of the CCdImages
    for (auto &ccdImage : ccdImageList) {
        const CcdImage &im = *ccdImage;
        auto visit = im.getVisit();
        auto chip = im.getCcdId();

        // check that the chip_indexed part was indeed assigned
        // (i.e. the reference visit was complete)
        if (_chipMap.find(chip) == _chipMap.end()) {
            LOGLS_WARN(_log, "Chip " << chip << " is missing in the reference exposure, expect troubles.");
            AstrometryTransformLinear norm = normalizeCoordinatesTransform(im.getImageFrame());
            _chipMap[chip] =
                    std::make_shared<SimplePolyMapping>(norm, AstrometryTransformPolynomial(chipOrder));
        }
        _mappings[ccdImage->getHashKey()] =
                std::make_unique<ChipVisitAstrometryMapping>(_chipMap[chip], _visitMap[visit]);
    }
    LOGLS_INFO(_log, "Got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                            << " visit mappings; holding chip " << constrainedChip << " fixed ("
                            << getTotalParameters() << " total parameters).");
    LOGLS_DEBUG(_log, "CcdImage map has " << _mappings.size() << " mappings, with "
                                          << _mappings.bucket_count() << " buckets and a load factor of "
                                          << _mappings.load_factor());
}

const AstrometryMapping *ConstrainedAstrometryModel::getMapping(CcdImage const &ccdImage) const {
    return findMapping(ccdImage);
}

/*! This routine decodes "DistortionsChip" and "DistortionsVisit" in
  whatToFit. If whatToFit contains "Distortions" and not
  Distortions<Something>, it is understood as both chips and
  visits. */
unsigned ConstrainedAstrometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    unsigned index = firstIndex;
    if (whatToFit.find("Distortions") == std::string::npos) {
        LOGLS_ERROR(_log, "assignIndices was called and Distortions is *not* in whatToFit");
        return 0;
    }
    // if we get here "Distortions" is in whatToFit
    _fittingChips = (whatToFit.find("DistortionsChip") != std::string::npos);
    _fittingVisits = (whatToFit.find("DistortionsVisit") != std::string::npos);
    // If nothing more than "Distortions" is specified, it means all:
    if ((!_fittingChips) && (!_fittingVisits)) {
        _fittingChips = _fittingVisits = true;
    }
    if (_fittingChips)
        for (auto &i : _chipMap) {
            i.second->setIndex(index);
            index += i.second->getNpar();
        }
    if (_fittingVisits)
        for (auto &i : _visitMap) {
            i.second->setIndex(index);
            index += i.second->getNpar();
        }
    // Tell the mappings which derivatives they will have to fill:
    for (auto &i : _mappings) {
        i.second->setWhatToFit(_fittingChips, _fittingVisits);
    }
    return index;
}

void ConstrainedAstrometryModel::offsetParams(Eigen::VectorXd const &delta) {
    if (_fittingChips)
        for (auto &i : _chipMap) {
            auto mapping = i.second.get();
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
    if (_fittingVisits)
        for (auto &i : _visitMap) {
            auto mapping = i.second.get();
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
}

void ConstrainedAstrometryModel::freezeErrorTransform() {
    for (auto i = _visitMap.begin(); i != _visitMap.end(); ++i) i->second->freezeErrorTransform();
    for (auto i = _chipMap.begin(); i != _chipMap.end(); ++i) i->second->freezeErrorTransform();
}

const AstrometryTransform &ConstrainedAstrometryModel::getChipTransform(CcdIdType const chip) const {
    auto chipp = _chipMap.find(chip);
    if (chipp == _chipMap.end()) {
        std::stringstream errMsg;
        errMsg << "No such chipId: " << chip << " among ";
        outputMapKeys(_chipMap, errMsg);
        std::cout << std::endl;
        throw pexExcept::InvalidParameterError(errMsg.str());
    }
    return chipp->second->getTransform();
}

// Array of visits involved in the solution.
std::vector<VisitIdType> ConstrainedAstrometryModel::getVisits() const {
    std::vector<VisitIdType> res;
    res.reserve(_visitMap.size());
    for (auto i = _visitMap.begin(); i != _visitMap.end(); ++i) res.push_back(i->first);
    return res;
}

const AstrometryTransform &ConstrainedAstrometryModel::getVisitTransform(VisitIdType const &visit) const {
    auto visitp = _visitMap.find(visit);
    if (visitp == _visitMap.end()) {
        std::stringstream errMsg;
        errMsg << "No such visitId: " << visit << " among ";
        outputMapKeys(_visitMap, errMsg);
        std::cout << std::endl;
        throw pexExcept::InvalidParameterError(errMsg.str());
    }
    return visitp->second->getTransform();
}

int ConstrainedAstrometryModel::getTotalParameters() const {
    int total = 0;
    for (auto &i : _chipMap) {
        total += i.second->getNpar();
    }
    for (auto &i : _visitMap) {
        total += i.second->getNpar();
    }
    return total;
}

std::shared_ptr<afw::geom::SkyWcs> ConstrainedAstrometryModel::makeSkyWcs(CcdImage const &ccdImage) const {
    auto proj = std::dynamic_pointer_cast<const TanRaDecToPixel>(getSkyToTangentPlane(ccdImage));
    jointcal::Point tangentPoint(proj->getTangentPoint());

    auto imageFrame = ccdImage.getImageFrame();
    auto pixelsToFocal = getChipTransform(ccdImage.getCcdId()).toAstMap(imageFrame);
    jointcal::Frame focalBox = getChipTransform(ccdImage.getCcdId()).apply(imageFrame, false);
    auto focalToIwc = getVisitTransform(ccdImage.getVisit()).toAstMap(focalBox);

    ast::Frame pixelFrame(2, "Domain=PIXELS");
    ast::Frame focalFrame(2, "Domain=FOCAL");
    ast::Frame iwcFrame(2, "Domain=IWC");

    // make a basic SkyWcs and extract the IWC portion
    auto iwcToSkyWcs = afw::geom::makeSkyWcs(
            afw::geom::Point2D(0, 0),
            afw::geom::SpherePoint(tangentPoint.x, tangentPoint.y, afw::geom::degrees),
            afw::geom::makeCdMatrix(1.0 * afw::geom::degrees, 0 * afw::geom::degrees, true));
    auto iwcToSkyMap = iwcToSkyWcs->getFrameDict()->getMapping("PIXELS", "SKY");
    auto skyFrame = iwcToSkyWcs->getFrameDict()->getFrame("SKY");

    ast::FrameDict frameDict(pixelFrame);
    frameDict.addFrame("PIXELS", *pixelsToFocal, focalFrame);
    frameDict.addFrame("FOCAL", *focalToIwc, iwcFrame);
    frameDict.addFrame("IWC", *iwcToSkyMap, *skyFrame);
    return std::make_shared<afw::geom::SkyWcs>(frameDict);
}

AstrometryMapping *ConstrainedAstrometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _mappings.find(ccdImage.getHashKey());
    if (i == _mappings.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "ConstrainedAstrometryModel cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
