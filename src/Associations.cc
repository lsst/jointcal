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

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/ListMatch.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/MeasuredStar.h"

#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"

#include "lsst/pex/exceptions.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/image/Calib.h"

#include "lsst/sphgeom/LonLat.h"
#include "lsst/sphgeom/Circle.h"
#include "lsst/sphgeom/ConvexPolygon.h"

namespace jointcal = lsst::jointcal;

namespace {
LOG_LOGGER _log = LOG_GET("lsst.jointcal.Associations");
}

namespace lsst {
namespace jointcal {

void Associations::createCcdImage(afw::table::SourceCatalog &catalog,
                                  const std::shared_ptr<lsst::afw::geom::SkyWcs>& wcs,
                                  const std::shared_ptr<lsst::afw::image::VisitInfo>& visitInfo,
                                  lsst::geom::Box2I const &bbox, std::string const &filter,
                                  const std::shared_ptr<afw::image::PhotoCalib>& photoCalib,
                                  const std::shared_ptr<afw::cameraGeom::Detector>& detector, int visit, int ccd,
                                  lsst::jointcal::JointcalControl const &control) {
    auto ccdImage = std::make_shared<CcdImage>(catalog, wcs, visitInfo, bbox, filter, photoCalib, detector,
                                               visit, ccd, control.sourceFluxField);
    ccdImageList.push_back(ccdImage);
}

void Associations::computeCommonTangentPoint() {
    std::vector<geom::SpherePoint> centers;
    centers.reserve(ccdImageList.size());
    for (auto const &ccdImage : ccdImageList) {
        centers.push_back(ccdImage->getBoresightRaDec());
    }
    auto commonTangentPoint = geom::averageSpherePoint(centers);
    LOGLS_DEBUG(_log, "Using common tangent point: " << commonTangentPoint.getPosition(geom::degrees));
    setCommonTangentPoint(commonTangentPoint.getPosition(geom::degrees));
}

void Associations::setCommonTangentPoint(lsst::geom::Point2D const &commonTangentPoint) {
    _commonTangentPoint = Point(commonTangentPoint.getX(), commonTangentPoint.getY());  // a jointcal::Point
    for (auto &ccdImage : ccdImageList) ccdImage->setCommonTangentPoint(_commonTangentPoint);
}

lsst::sphgeom::Circle Associations::computeBoundingCircle() const {
    // Compute the frame on the common tangent plane that contains all input images.
    Frame tangentPlaneFrame;

    for (auto const &ccdImage : ccdImageList) {
        Frame CTPFrame = ccdImage->getPixelToCommonTangentPlane()->apply(ccdImage->getImageFrame(), false);
        if (tangentPlaneFrame.getArea() == 0)
            tangentPlaneFrame = CTPFrame;
        else
            tangentPlaneFrame += CTPFrame;
    }

    // Convert tangent plane coordinates to RaDec.
    AstrometryTransformLinear identity;
    TanPixelToRaDec commonTangentPlaneToRaDec(identity, _commonTangentPoint);
    Frame raDecFrame = commonTangentPlaneToRaDec.apply(tangentPlaneFrame, false);

    std::vector<sphgeom::UnitVector3d> points;
    points.reserve(4);
    // raDecFrame is in on-sky (RA,Dec) degrees stored as an x/y box:
    // the on-sky bounding box it represents is given by the corners of that box.
    LOGLS_INFO(_log, "Computed tangent plane box for data (degrees): " << raDecFrame);
    points.emplace_back(sphgeom::LonLat::fromDegrees(raDecFrame.xMin, raDecFrame.yMin));
    points.emplace_back(sphgeom::LonLat::fromDegrees(raDecFrame.xMax, raDecFrame.yMin));
    points.emplace_back(sphgeom::LonLat::fromDegrees(raDecFrame.xMin, raDecFrame.yMax));
    points.emplace_back(sphgeom::LonLat::fromDegrees(raDecFrame.xMax, raDecFrame.yMax));

    return sphgeom::ConvexPolygon::convexHull(points).getBoundingCircle();
}

void Associations::associateCatalogs(const double matchCutInArcSec, const bool useFittedList,
                                     const bool enlargeFittedList) {
    // clear reference stars
    refStarList.clear();

    // clear measurement counts and associations to refstars, but keep fittedStars themselves.
    for (auto &item : fittedStarList) {
        item->clearBeforeAssoc();
    }
    // clear fitted stars
    if (!useFittedList) fittedStarList.clear();

    for (auto &ccdImage : ccdImageList) {
        std::shared_ptr<AstrometryTransform> toCommonTangentPlane = ccdImage->getPixelToCommonTangentPlane();

        // Clear the catalog to fit and copy the whole catalog into it.
        // This allows reassociating from scratch after a fit.
        ccdImage->resetCatalogForFit();
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();

        // Associate with previous lists.
        /* To speed up the match (more precisely the contruction of the FastFinder), select in the
         fittedStarList the objects that are within reach of the current ccdImage */
        Frame ccdImageFrameCPT = toCommonTangentPlane->apply(ccdImage->getImageFrame(), false);
        ccdImageFrameCPT = ccdImageFrameCPT.rescale(1.10);  // add 10 % margin.
        // We cannot use FittedStarList::ExtractInFrame, because it does an actual copy, which we don't want
        // here: we want the pointers in the StarMatch to refer to fittedStarList elements.
        FittedStarList toMatch;

        for (auto const &fittedStar : fittedStarList) {
            if (ccdImageFrameCPT.inFrame(*fittedStar)) {
                toMatch.push_back(fittedStar);
            }
        }

        // divide by 3600 because coordinates in CTP are in degrees.
        auto starMatchList = listMatchCollect(Measured2Base(catalog), Fitted2Base(toMatch),
                                              toCommonTangentPlane.get(), matchCutInArcSec / 3600.);

        /* should check what this removeAmbiguities does... */
        LOGLS_DEBUG(_log, "Measured-to-Fitted matches before removing ambiguities " << starMatchList->size());
        starMatchList->removeAmbiguities(*toCommonTangentPlane);
        LOGLS_DEBUG(_log, "Measured-to-Fitted matches after removing ambiguities " << starMatchList->size());

        // Associate MeasuredStar -> FittedStar using the surviving matches.

        int matchedCount = 0;
        for (auto const &starMatch : *starMatchList) {
            auto bs = starMatch.s1;
            auto ms_const = std::dynamic_pointer_cast<const MeasuredStar>(bs);
            auto ms = std::const_pointer_cast<MeasuredStar>(ms_const);
            auto bs2 = starMatch.s2;
            auto fs_const = std::dynamic_pointer_cast<const FittedStar>(bs2);
            auto fs = std::const_pointer_cast<FittedStar>(fs_const);
            ms->setFittedStar(fs);
            matchedCount++;
        }
        LOGLS_DEBUG(_log, "Matched " << matchedCount << " objects in " << ccdImage->getName());

        // add unmatched objets to FittedStarList
        int unMatchedCount = 0;
        for (auto const &mstar : catalog) {
            // to check if it was matched, just check if it has a fittedStar Pointer assigned
            if (mstar->getFittedStar()) continue;
            if (enlargeFittedList) {
                auto fs = std::make_shared<FittedStar>(*mstar);
                // transform coordinates to CommonTangentPlane
                toCommonTangentPlane->transformPosAndErrors(*fs, *fs);
                fittedStarList.push_back(fs);
                mstar->setFittedStar(fs);
            }
            unMatchedCount++;
        }
        LOGLS_DEBUG(_log, "Unmatched objects: " << unMatchedCount);
    }  // end of loop on CcdImages

    // !!!!!!!!!!!!!!!!!
    // TODO: DO WE REALLY NEED THIS???
    // Why do we need to do this, instead of directly computing them in normalizeFittedStars?
    // What makes the magnitudes special here?
    // !!!!!!!!!!!!!!!!!
    // assignMags();
}

void Associations::collectRefStars(afw::table::SimpleCatalog &refCat, geom::Angle matchCut,
                                   std::string const &fluxField, float refCoordinateErr,
                                   bool rejectBadFluxes) {
    if (refCat.size() == 0) {
        throw(LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " reference catalog is empty : stop here "));
    }

    afw::table::CoordKey coordKey = refCat.getSchema()["coord"];
    // Handle reference catalogs that don't have position errors.
    afw::table::Key<float> raErrKey;
    afw::table::Key<float> decErrKey;
    if (std::isnan(refCoordinateErr)) {
        raErrKey = refCat.getSchema()["coord_raErr"];
        decErrKey = refCat.getSchema()["coord_decErr"];
    }

    auto fluxKey = refCat.getSchema().find<double>(fluxField).key;
    // Handle reference catalogs that don't have flux errors.
    afw::table::Key<double> fluxErrKey;
    try {
        fluxErrKey = refCat.getSchema().find<double>(fluxField + "Err").key;
    } catch (pex::exceptions::NotFoundError &) {
        LOGLS_WARN(_log, "Flux error field ("
                                 << fluxField << "Err"
                                 << ") not found in reference catalog. Not using ref flux errors.");
    }

    // Handle reference catalogs that don't have proper motion & error
    afw::table::Key<geom::Angle> pmRaKey, pmDecKey;
    afw::table::Key<float> pmRaErrKey, pmDecErrKey, pmRaDecCovKey;
    try {
        pmRaKey = refCat.getSchema().find<geom::Angle>("pm_ra").key;
        pmDecKey = refCat.getSchema().find<geom::Angle>("pm_dec").key;
    } catch (pex::exceptions::NotFoundError &ex) {
        LOGLS_WARN(_log, "Not loading proper motions: (pm_ra,pm_dec) fields not found in reference catalog.");
    } catch (pex::exceptions::TypeError &ex) {
        LOGLS_WARN(_log, "Not loading proper motions: RA/Dec proper motion values must be `geom:Angle`: "
                                 << ex.what());
    }
    try {
        pmRaErrKey = refCat.getSchema().find<float>("pm_raErr").key;
        pmDecErrKey = refCat.getSchema().find<float>("pm_decErr").key;
    } catch (pex::exceptions::NotFoundError &ex) {
        LOGLS_WARN(_log, "Not loading proper motions: error fields not available: " << ex.what());
    }

    // TODO: we aren't getting covariances from Gaia yet, so maybe ignore this for now?
    try {
        pmRaDecCovKey = refCat.getSchema().find<float>("pm_ra_Dec_Cov").key;
    } catch (pex::exceptions::NotFoundError &ex) {
        LOGLS_WARN(_log, "No ra/dec proper motion covariances in refcat: " << ex.what());
    }

    refStarList.clear();
    for (size_t i = 0; i < refCat.size(); i++) {
        auto const &record = refCat.get(i);

        auto coord = record->get(coordKey);
        double flux = record->get(fluxKey);
        double fluxErr;
        if (fluxErrKey.isValid()) {
            fluxErr = record->get(fluxErrKey);
        } else {
            fluxErr = std::numeric_limits<double>::quiet_NaN();
        }
        double ra = lsst::geom::radToDeg(coord.getLongitude());
        double dec = lsst::geom::radToDeg(coord.getLatitude());
        auto star = std::make_shared<RefStar>(ra, dec, flux, fluxErr);

        if (std::isnan(refCoordinateErr)) {
            // refcat errors are unitless but stored as radians: convert to deg**2
            star->vx = std::pow(lsst::geom::radToDeg(record->get(raErrKey)), 2);
            star->vy = std::pow(lsst::geom::radToDeg(record->get(decErrKey)), 2);
        } else {
            // Convert the fake errors from mas to deg**2
            star->vx = std::pow(refCoordinateErr / 1000. / 3600. / std::cos(coord.getLatitude()), 2);
            star->vy = std::pow(refCoordinateErr / 1000. / 3600., 2);
        }

        if (pmRaKey.isValid()) {
            if (pmRaDecCovKey.isValid()) {
                star->setProperMotion(std::make_unique<ProperMotion const>(
                        record->get(pmRaKey).asRadians(), record->get(pmDecKey).asRadians(),
                        record->get(pmRaErrKey), record->get(pmDecErrKey), record->get(pmRaDecCovKey)));
            } else {
                star->setProperMotion(std::make_unique<ProperMotion const>(
                        record->get(pmRaKey).asRadians(), record->get(pmDecKey).asRadians(),
                        record->get(pmRaErrKey), record->get(pmDecErrKey)));
            }
        }

        // TODO: cook up a covariance as none of our current refcats have it
        star->vxy = 0.;

        // Reject sources with non-finite fluxes and flux errors, and fluxErr=0 (which gives chi2=inf).
        if (rejectBadFluxes && (!std::isfinite(flux) || !std::isfinite(fluxErr) || fluxErr <= 0)) continue;
        refStarList.push_back(star);
    }

    // project on CTP (i.e. RaDec2CTP), in degrees
    AstrometryTransformLinear identity;
    TanRaDecToPixel raDecToCommonTangentPlane(identity, _commonTangentPoint);

    associateRefStars(matchCut.asArcseconds(), &raDecToCommonTangentPlane);
}

void Associations::associateRefStars(double matchCutInArcSec, const AstrometryTransform *transform) {
    // associate with FittedStars
    // 3600 because coordinates are in degrees (in CTP).
    auto starMatchList = listMatchCollect(Ref2Base(refStarList), Fitted2Base(fittedStarList), transform,
                                          matchCutInArcSec / 3600.);

    LOGLS_DEBUG(_log, "Refcat matches before removing ambiguities " << starMatchList->size());
    starMatchList->removeAmbiguities(*transform);
    LOGLS_DEBUG(_log, "Refcat matches after removing ambiguities " << starMatchList->size());

    // actually associate things
    for (auto const &starMatch : *starMatchList) {
        const BaseStar &bs = *starMatch.s1;
        const auto &rs_const = dynamic_cast<const RefStar &>(bs);
        auto &rs = const_cast<RefStar &>(rs_const);
        const BaseStar &bs2 = *starMatch.s2;
        const auto &fs_const = dynamic_cast<const FittedStar &>(bs2);
        auto &fs = const_cast<FittedStar &>(fs_const);
        // rs->setFittedStar(*fs);
        fs.setRefStar(&rs);
    }

    LOGLS_INFO(_log,
               "Associated " << starMatchList->size() << " reference stars among " << refStarList.size());
}

void Associations::prepareFittedStars(int minMeasurements) {
    selectFittedStars(minMeasurements);
    normalizeFittedStars();
}

void Associations::selectFittedStars(int minMeasurements) {
    LOGLS_INFO(_log, "Fitted stars before measurement # cut: " << fittedStarList.size());

    int totalMeasured = 0, validMeasured = 0;

    // first pass: remove objects that have less than a certain number of measurements.
    for (auto const &ccdImage : ccdImageList) {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        // Iteration happens internal to the loop, as we may delete measuredStars from catalog.
        for (auto mi = catalog.begin(); mi != catalog.end();) {
            MeasuredStar &mstar = **mi;
            ++totalMeasured;

            auto fittedStar = mstar.getFittedStar();
            // measuredStar has no fittedStar: move on.
            if (fittedStar == nullptr) {
                ++mi;
                continue;
            }

            // keep FittedStars which either have a minimum number of
            // measurements, or are matched to a RefStar
            if (!fittedStar->getRefStar() && fittedStar->getMeasurementCount() < minMeasurements) {
                fittedStar->getMeasurementCount()--;
                mi = catalog.erase(mi);  // mi now points to the next measuredStar.
            } else {
                ++validMeasured;
                ++mi;
            }
        }  // end loop on objects in catalog
    }      // end loop on catalogs

    // now FittedStars with less than minMeasurements should have zero measurementCount.
    for (auto fi = fittedStarList.begin(); fi != fittedStarList.end();) {
        if ((*fi)->getMeasurementCount() == 0) {
            fi = fittedStarList.erase(fi);
        } else {
            ++fi;
        }
    }

    LOGLS_INFO(_log, "Fitted stars after measurement # cut: " << fittedStarList.size());
    LOGLS_INFO(_log, "Total, valid number of Measured stars: " << totalMeasured << ", " << validMeasured);
}

void Associations::normalizeFittedStars() {
    // Clear positions in order to take the average of the measuredStars.
    for (auto &fittedStar : fittedStarList) {
        fittedStar->x = 0.0;
        fittedStar->y = 0.0;
        fittedStar->setFlux(0.0);
        fittedStar->getMag() = 0.0;
    }

    // Iterate over measuredStars to add their values into their fittedStars
    for (auto const &ccdImage : ccdImageList) {
        std::shared_ptr<AstrometryTransform> toCommonTangentPlane = ccdImage->getPixelToCommonTangentPlane();
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (auto &mi : catalog) {
            auto fittedStar = mi->getFittedStar();
            if (fittedStar == nullptr)
                throw(LSST_EXCEPT(
                        pex::exceptions::RuntimeError,
                        "All measuredStars must have a fittedStar: did you call selectFittedStars()?"));
            auto point = toCommonTangentPlane->apply(*mi);
            fittedStar->x += point.x;
            fittedStar->y += point.y;
            fittedStar->getFlux() += mi->getFlux();
        }
    }

    _maxMeasuredStars = 0;
    for (auto &fi : fittedStarList) {
        auto measurementCount = fi->getMeasurementCount();
        _maxMeasuredStars += measurementCount;
        fi->x /= measurementCount;
        fi->y /= measurementCount;
        fi->getFlux() /= measurementCount;
        fi->getMag() = utils::nanojanskyToABMagnitude(fi->getFlux());
    }
}

void Associations::cleanFittedStars() {
    auto iter = fittedStarList.begin();
    auto end = fittedStarList.end();
    size_t count = 0;
    while (iter != end) {
        auto fittedStar = *iter;
        if (fittedStar->getMeasurementCount() == 0) {
            LOGLS_TRACE(_log, "Deleting FittedStar (has no measuredStars): " << *fittedStar);
            iter = fittedStarList.erase(iter);
            count++;
        } else {
            ++iter;
        }
    }
    if (count > 0) {
        LOGLS_INFO(_log, "Removed " << count << " fittedStars that had no associated measuredStar.");
    }
}

void Associations::assignMags() {
    for (auto const &ccdImage : ccdImageList) {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (auto const &mstar : catalog) {
            auto fstar = mstar->getFittedStar();
            if (!fstar) continue;
            fstar->addMagMeasurement(mstar->getMag(), mstar->getMagWeight());
        }
    }
}

void Associations::deprojectFittedStars() {
    // By default, Associations::fittedStarList is expressed on the Associations::commonTangentPlane.
    // For AstrometryFit, we need it on the sky.
    if (!fittedStarList.inTangentPlaneCoordinates) {
        LOGLS_WARN(_log,
                   "DeprojectFittedStars: Fitted stars are already in sidereal coordinates, nothing done ");
        return;
    }

    TanPixelToRaDec ctp2Sky(AstrometryTransformLinear(), getCommonTangentPoint());
    fittedStarList.applyTransform(ctp2Sky);
    fittedStarList.inTangentPlaneCoordinates = false;
}

int Associations::nCcdImagesValidForFit() const {
    return std::count_if(ccdImageList.begin(), ccdImageList.end(), [](std::shared_ptr<CcdImage> const &item) {
        return item->getCatalogForFit().size() > 0;
    });
}

size_t Associations::nFittedStarsWithAssociatedRefStar() const {
    size_t count = 0;
    for (auto const &fittedStar : fittedStarList) {
        if ((fittedStar != nullptr) & (fittedStar->getRefStar() != nullptr)) count++;
    }
    return count;
}

#ifdef TODO
void Associations::collectMCStars(int realization) {
    CcdImageIterator I;
    StarMatchIterator smI;

    for (I = ccdImageList.begin(); I != ccdImageList.end(); I++) {
        CcdImage &ccdImage = **I;
        string dbimdir = ccdImage.Dir();
        string mctruth = dbimdir + "/mc/mctruth.list";

        if (realization >= 0) {
            stringstream sstrm;
            sstrm << dbimdir << "/mc/mctruth_" << realization << ".list";
            mctruth = sstrm.str();
        }

        AstrometryTransformIdentity gti;
        MeasuredStarList &catalog = ccdImage.getCatalogForFit();

        //      BaseStarWithErrorList mctruthlist(mctruth);
        DicStarList mctruthlist(mctruth);
        auto starMatchList =
                listMatchCollect(Measured2Base(catalog), Dic2Base(mctruthlist), &gti, 1. /* pixel ? */);
        if (starMatchList)
            for (smI = starMatchList->begin(); smI != starMatchList->end(); smI++) {
                StarMatch &sm = *smI;
                BaseStar *bs = sm.s1;
                MeasuredStar *mstar = dynamic_cast<MeasuredStar *>(bs);
                bs = sm.s2;
                DicStar *dstar = dynamic_cast<DicStar *>(bs);
                std::unique_ptr<BaseStarWithError> mcstar(new BaseStarWithError(*bs));
                mcstar->GetMCInfo().iflux = dstar->getval("iflux");
                mcstar->GetMCInfo().tflux = dstar->getval("sflux");
                /*
                mstar->SetMCTruth(mcstar);
                mstar->SetMCMeas(mcstar);
                */
            }
        else
            LOGLS_FATAL(_log, "CollectMCStars Unable to match MCTruth w/ catalog!");
    }
}

#endif /* TODO */
}  // namespace jointcal
}  // namespace lsst
