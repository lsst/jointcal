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

#include <memory>
#include <string>

#include "astshim.h"
#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimpleAstrometryModel.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/AstrometryTransform.h"

namespace lsst {
namespace jointcal {

SimpleAstrometryModel::SimpleAstrometryModel(CcdImageList const &ccdImageList,
                                             const std::shared_ptr<ProjectionHandler const> projectionHandler,
                                             bool initFromWcs, unsigned nNotFit, unsigned order)
        : AstrometryModel(LOG_GET("jointcal.SimpleAstrometryModel")),
          _skyToTangentPlane(projectionHandler)

{
    unsigned count = 0;

    for (auto i = ccdImageList.cbegin(); i != ccdImageList.cend(); ++i, ++count) {
        const CcdImage &im = **i;
        if (count < nNotFit) {
            std::unique_ptr<SimpleAstrometryMapping> id(
                    new SimpleAstrometryMapping(AstrometryTransformIdentity()));
            id->setIndex(-1);  // non sense, because it has no parameters
            _myMap[im.getHashKey()] = std::move(id);
        } else
        // Given how AssignIndices works, only the SimplePolyMappings
        // will actually be fitted, as nNotFit requests.
        {
            // first check that there are enough measurements for the requested polynomial order.
            size_t nObj = im.getCatalogForFit().size();
            if (nObj == 0) {
                LOGLS_WARN(_log, "Empty catalog from image: " << im.getName());
                continue;
            }
            AstrometryTransformPolynomial pol(order);
            if (pol.getOrder() > 0)  // if not, it cannot be decreased
            {
                while (unsigned(pol.getNpar()) > 2 * nObj) {
                    LOGLS_WARN(_log, "Reducing polynomial order from "
                                             << pol.getOrder() << ", due to too few sources (" << nObj
                                             << " vs. " << pol.getNpar() << " parameters)");
                    pol.setOrder(pol.getOrder() - 1);
                }
            }
            /* We have to center and normalize the coordinates so that
               the fit matrix is not too ill-conditionned. Basically, x
               and y in pixels are mapped to [-1,1]. When the
               transformation of SimplePolyMapping transformation is
               accessed, the combination of the normalization and the
               fitted transformation is returned, so that the trick
               remains hidden
             */
            const Frame &frame = im.getImageFrame();
            AstrometryTransformLinear shiftAndNormalize = normalizeCoordinatesTransform(frame);
            if (initFromWcs) {
                pol = AstrometryTransformPolynomial(im.getPixelToTangentPlane().get(), frame, order);
                pol = pol * shiftAndNormalize.inverted();
            }
            _myMap[im.getHashKey()] =
                    std::unique_ptr<SimpleAstrometryMapping>(new SimplePolyMapping(shiftAndNormalize, pol));
        }
    }
}

const AstrometryMapping *SimpleAstrometryModel::getMapping(CcdImage const &ccdImage) const {
    return findMapping(ccdImage);
}

unsigned SimpleAstrometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    if (whatToFit.find("Distortions") == std::string::npos) {
        LOGLS_ERROR(_log, "AssignIndices was called and Distortions is *not* in whatToFit.");
        return 0;
    }
    unsigned index = firstIndex;
    for (auto i = _myMap.begin(); i != _myMap.end(); ++i) {
        SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
        if (!p) continue;  // it should be AstrometryTransformIdentity
        p->setIndex(index);
        index += p->getNpar();
    }
    return index;
}

void SimpleAstrometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

void SimpleAstrometryModel::freezeErrorTransform() {
    for (auto &i : _myMap) i.second->freezeErrorTransform();
}

int SimpleAstrometryModel::getTotalParameters() const {
    int total = 0;
    for (auto &i : _myMap) {
        total += i.second->getNpar();
    }
    return total;
}

const AstrometryTransform &SimpleAstrometryModel::getTransform(CcdImage const &ccdImage) const {
    return dynamic_cast<const SimplePolyMapping *>(findMapping(ccdImage))->getTransform();
}

std::shared_ptr<afw::geom::SkyWcs> SimpleAstrometryModel::makeSkyWcs(CcdImage const &ccdImage) const {
    auto proj = std::dynamic_pointer_cast<const TanRaDecToPixel>(getSkyToTangentPlane(ccdImage));
    jointcal::Point tangentPoint(proj->getTangentPoint());

    auto polyMap = getTransform(ccdImage).toAstMap(ccdImage.getImageFrame());
    ast::Frame pixelFrame(2, "Domain=PIXELS");
    ast::Frame iwcFrame(2, "Domain=IWC");

    // make a basic SkyWcs and extract the IWC portion
    auto iwcToSkyWcs = afw::geom::makeSkyWcs(
            afw::geom::Point2D(0, 0),
            afw::geom::SpherePoint(tangentPoint.x, tangentPoint.y, afw::geom::degrees),
            afw::geom::makeCdMatrix(1.0 * afw::geom::degrees, 0 * afw::geom::degrees, true));
    auto iwcToSkyMap = iwcToSkyWcs->getFrameDict()->getMapping("PIXELS", "SKY");
    auto skyFrame = iwcToSkyWcs->getFrameDict()->getFrame("SKY");

    ast::FrameDict frameDict(pixelFrame, *polyMap, iwcFrame);
    frameDict.addFrame("IWC", *iwcToSkyMap, *skyFrame);
    return std::make_shared<afw::geom::SkyWcs>(frameDict);
}

AstrometryMapping *SimpleAstrometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimpleAstrometryModel cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
