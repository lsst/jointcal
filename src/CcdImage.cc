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

#include <assert.h>
#include <string>
#include <sstream>
#include <cmath>

#include "lsst/afw/cameraGeom/CameraSys.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Point.h"

#include "lsst/log/Log.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Point.h"

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.CcdImage");
}

namespace lsst {
namespace jointcal {

std::ostream &operator<<(std::ostream &out, CcdImageKey const &key) {
    out << "(visit: " << key.visit << ", detector: " << key.ccd << ")";
    return out;
}

void CcdImage::loadCatalog(afw::table::SourceCatalog const &catalog, std::string const &fluxField) {
    auto xKey = catalog.getSchema().find<double>("slot_Centroid_x").key;
    auto yKey = catalog.getSchema().find<double>("slot_Centroid_y").key;
    auto xsKey = catalog.getSchema().find<float>("slot_Centroid_xErr").key;
    auto ysKey = catalog.getSchema().find<float>("slot_Centroid_yErr").key;
    auto mxxKey = catalog.getSchema().find<double>("slot_Shape_xx").key;
    auto myyKey = catalog.getSchema().find<double>("slot_Shape_yy").key;
    auto mxyKey = catalog.getSchema().find<double>("slot_Shape_xy").key;
    auto instFluxKey = catalog.getSchema().find<double>(fluxField + "_instFlux").key;
    auto instFluxErrKey = catalog.getSchema().find<double>(fluxField + "_instFluxErr").key;

    auto transform = _detector->getTransform(afw::cameraGeom::PIXELS, afw::cameraGeom::FOCAL_PLANE);

    _wholeCatalog.clear();
    for (auto const &record : catalog) {
        auto ms = std::make_shared<MeasuredStar>();
        ms->setId(record.getId());
        ms->x = record.get(xKey);
        ms->y = record.get(yKey);
        ms->vx = std::pow(record.get(xsKey), 2);
        ms->vy = std::pow(record.get(ysKey), 2);
        auto pointFocal = transform->applyForward(record.getCentroid());
        ms->setXFocal(pointFocal.getX());
        ms->setYFocal(pointFocal.getY());
        /* the xy covariance is not provided in the input catalog: we
        cook it up from the x and y position variance and the shape
         measurements: */
        double mxx = record.get(mxxKey);
        double myy = record.get(myyKey);
        double mxy = record.get(mxyKey);
        ms->vxy = mxy * (ms->vx + ms->vy) / (mxx + myy);
        if (std::isnan(ms->vxy) || ms->vx < 0 || ms->vy < 0 || (ms->vxy * ms->vxy) > (ms->vx * ms->vy)) {
            LOGLS_WARN(_log, "Bad source detected during loadCatalog id: "
                                     << ms->getId() << " with vx,vy: " << ms->vx << "," << ms->vy
                                     << " vxy^2: " << ms->vxy * ms->vxy << " vx*vy: " << ms->vx * ms->vy);
            continue;
        }
        ms->setInstFluxAndErr(record.get(instFluxKey), record.get(instFluxErrKey));
        // TODO: the below lines will be less clumsy once DM-4044 is cleaned up and we can say:
        // TODO: instFluxToNanojansky(ms->getInstFlux(), ms) (because ms will be derived from
        // geom::Point).
        geom::Point<double, 2> point(ms->x, ms->y);
        auto flux = _photoCalib->instFluxToNanojansky(ms->getInstFlux(), ms->getInstFluxErr(), point);
        ms->setFlux(flux.value);
        ms->setFluxErr(flux.error);
        auto mag = _photoCalib->instFluxToMagnitude(ms->getInstFlux(), ms->getInstFluxErr(), point);
        ms->getMag() = mag.value;
        ms->setMagErr(mag.error);
        ms->setCcdImage(this);
        _wholeCatalog.push_back(std::move(ms));
    }
    _wholeCatalog.setCcdImage(this);
}

CcdImage::CcdImage(afw::table::SourceCatalog &catalog, std::shared_ptr<lsst::afw::geom::SkyWcs> wcs,
                   std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, geom::Box2I const &bbox,
                   std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib,
                   std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccdId,
                   std::string const &fluxField)
        : _ccdId(ccdId), _visit(visit), _photoCalib(photoCalib), _detector(detector), _filter(filter) {
    loadCatalog(catalog, fluxField);

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    _imageFrame = Frame(lowerLeft, upperRight);

    _readWcs = std::make_shared<AstrometryTransformSkyWcs>(wcs);

    std::stringstream out;
    out << visit << "_" << ccdId;
    _name = out.str();

    _boresightRaDec = visitInfo->getBoresightRaDec();
    _airMass = visitInfo->getBoresightAirmass();
    _epoch = visitInfo->getDate().get(lsst::daf::base::DateTime::EPOCH);
    double latitude = visitInfo->getObservatory().getLatitude();
    _lstObs = visitInfo->getEra();
    _hourAngle = visitInfo->getBoresightHourAngle();

    // Some cameras don't manage ERA (and thus Hour Angle) properly, so it's going to be NaN.
    // Because we need the refraction vector later, go with 0 HA to prevent crashes on that NaN.
    if (std::isnan(_hourAngle) == true) {
        _hourAngle = 0;
    }

    if (_airMass == 1)
        _sinEta = _cosEta = _tanZ = 0;
    else {
        double cosz = 1. / _airMass;
        double sinz = std::sqrt(1 - cosz * cosz);  // astronomers usually observe above the horizon
        _tanZ = sinz / cosz;
        // TODO: as part of DM-12473, we can remove all of this and just call _visitInfo.getParallacticAngle()
        double dec = _boresightRaDec.getLatitude();
        // x/y components of refraction angle, eta.]
        double yEta = std::sin(_hourAngle);
        double xEta = std::cos(dec) * std::tan(latitude) - std::sin(dec) * std::cos(_hourAngle);
        double eta = std::atan2(yEta, xEta);
        _sinEta = std::sin(eta);
        _cosEta = std::cos(eta);
    }
}

std::pair<int, int> CcdImage::countStars() const {
    int measuredStars = 0;
    int refStars = 0;
    for (auto const &measuredStar : _catalogForFit) {
        if (measuredStar->isValid()) {
            measuredStars++;
        }
        if ((measuredStar->getFittedStar() != nullptr) &&
            (measuredStar->getFittedStar()->getRefStar() != nullptr)) {
            refStars++;
        }
    }
    return std::make_pair(measuredStars, refStars);
}

void CcdImage::setCommonTangentPoint(Point const &commonTangentPoint) {
    _commonTangentPoint = commonTangentPoint;

    auto const crval = _readWcs->getSkyWcs()->getSkyOrigin();
    jointcal::Point tangentPoint(crval[0].asDegrees(), crval[1].asDegrees());

    /* we don't assume here that we know the internals of TanPixelToRaDec:
       to construct pix->TP, we do pix->sky->TP, although pix->sky
       actually goes through TP */
    AstrometryTransformLinear identity;
    TanRaDecToPixel raDecToTangentPlane(identity, tangentPoint);
    _pixelToTangentPlane = compose(raDecToTangentPlane, *_readWcs);
    TanPixelToRaDec CommonTangentPlane2RaDec(identity, commonTangentPoint);
    _commonTangentPlaneToTangentPlane = compose(raDecToTangentPlane, CommonTangentPlane2RaDec);

    // jump from one TP to an other:
    TanRaDecToPixel raDecToCommonTangentPlane(identity, commonTangentPoint);
    TanPixelToRaDec TangentPlaneToRaDec(identity, tangentPoint);
    _tangentPlaneToCommonTangentPlane = compose(raDecToCommonTangentPlane, TangentPlaneToRaDec);
    _skyToTangentPlane.reset(new TanRaDecToPixel(identity, tangentPoint));

    // this one is needed for matches :
    _pixelToCommonTangentPlane = compose(raDecToCommonTangentPlane, *_readWcs);
}
}  // namespace jointcal
}  // namespace lsst
