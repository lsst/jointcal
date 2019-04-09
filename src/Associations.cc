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
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image/Calib.h"

namespace jointcal = lsst::jointcal;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.Associations");
}

namespace lsst {
namespace jointcal {

void Associations::createCcdImage(afw::table::SourceCatalog &catalog,
                                  std::shared_ptr<lsst::afw::geom::SkyWcs> wcs,
                                  std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo,
                                  lsst::afw::geom::Box2I const &bbox, std::string const &filter,
                                  std::shared_ptr<afw::image::PhotoCalib> photoCalib,
                                  std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
                                  lsst::jointcal::JointcalControl const &control) {
    auto ccdImage = std::make_shared<CcdImage>(catalog, wcs, visitInfo, bbox, filter, photoCalib, detector,
                                               visit, ccd, control.sourceFluxField);
    ccdImageList.push_back(ccdImage);
    LOGLS_DEBUG(_log, "Catalog " << ccdImage->getName() << " has " << ccdImage->getWholeCatalog().size()
                                 << " objects.");
}

void Associations::computeCommonTangentPoint() {
    std::vector<afw::geom::SpherePoint> centers;
    centers.reserve(ccdImageList.size());
    for (auto const &ccdImage : ccdImageList) {
        centers.push_back(ccdImage->getBoresightRaDec());
    }
    auto commonTangentPoint = afw::geom::averageSpherePoint(centers);
    LOGLS_DEBUG(_log, "Using common tangent point: " << commonTangentPoint.getPosition(afw::geom::degrees));
    setCommonTangentPoint(commonTangentPoint.getPosition(afw::geom::degrees));
}

void Associations::setCommonTangentPoint(lsst::afw::geom::Point2D const &commonTangentPoint) {
    _commonTangentPoint = Point(commonTangentPoint.getX(), commonTangentPoint.getY());  // a jointcal::Point
    for (auto &ccdImage : ccdImageList) ccdImage->setCommonTangentPoint(_commonTangentPoint);
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
        LOGLS_INFO(_log, "Matched " << matchedCount << " objects in " << ccdImage->getName());

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
        LOGLS_INFO(_log, "Unmatched objects: " << unMatchedCount);
    }  // end of loop on CcdImages

    // !!!!!!!!!!!!!!!!!
    // TODO: DO WE REALLY NEED THIS???
    // Why do we need to do this, instead of directly computing them in normalizeFittedStars?
    // What makes the magnitudes special here?
    // !!!!!!!!!!!!!!!!!
    // assignMags();
}

void Associations::collectRefStars(afw::table::SimpleCatalog &refCat, afw::geom::Angle matchCut,
                                   std::string const &fluxField, bool rejectBadFluxes) {
    if (refCat.size() == 0) {
        throw(LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " reference catalog is empty : stop here "));
    }

    afw::table::CoordKey coordKey = refCat.getSchema()["coord"];
    auto fluxKey = refCat.getSchema().find<double>(fluxField).key;
    // Don't blow up if the reference catalog doesn't contain errors.
    afw::table::Key<double> fluxErrKey;
    try {
        fluxErrKey = refCat.getSchema().find<double>(fluxField + "Err").key;
    } catch (pex::exceptions::NotFoundError &) {
        LOGLS_WARN(_log, "Flux error field ("
                                 << fluxField << "Err"
                                 << ") not found in reference catalog. Not using ref flux errors.");
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
        double ra = lsst::afw::geom::radToDeg(coord.getLongitude());
        double dec = lsst::afw::geom::radToDeg(coord.getLatitude());
        auto star = std::make_shared<RefStar>(ra, dec, flux, fluxErr);

        // TODO DM-10826: RefCats aren't guaranteed to have position errors.
        // TODO: Need to devise a way to check whether the refCat has position errors
        // TODO: and use them instead, if available.
        // cook up errors: 100 mas per cooordinate
        star->vx = std::pow(0.1 / 3600 / cos(coord.getLatitude()), 2);
        star->vy = std::pow(0.1 / 3600, 2);
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

const lsst::afw::geom::Box2D Associations::getRaDecBBox() {
    // compute the frame on the CTP that contains all input images
    Frame tangentPlaneFrame;

    for (auto const &ccdImage : ccdImageList) {
        Frame CTPFrame = ccdImage->getPixelToCommonTangentPlane()->apply(ccdImage->getImageFrame(), false);
        if (tangentPlaneFrame.getArea() == 0)
            tangentPlaneFrame = CTPFrame;
        else
            tangentPlaneFrame += CTPFrame;
    }

    // convert tangent plane coordinates to RaDec:
    AstrometryTransformLinear identity;
    TanPixelToRaDec commonTangentPlaneToRaDec(identity, _commonTangentPoint);
    Frame raDecFrame = commonTangentPlaneToRaDec.apply(tangentPlaneFrame, false);

    lsst::afw::geom::Point<double> min(raDecFrame.xMin, raDecFrame.yMin);
    lsst::afw::geom::Point<double> max(raDecFrame.xMax, raDecFrame.yMax);
    lsst::afw::geom::Box2D box(min, max);

    return box;
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
        const RefStar &rs_const = dynamic_cast<const RefStar &>(bs);
        RefStar &rs = const_cast<RefStar &>(rs_const);
        const BaseStar &bs2 = *starMatch.s2;
        const FittedStar &fs_const = dynamic_cast<const FittedStar &>(bs2);
        FittedStar &fs = const_cast<FittedStar &>(fs_const);
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

    std::size_t totalMeasured = 0, validMeasured = 0;

    // first pass: remove objects that have less than a certain number of measurements.
    for (auto const &ccdImage : ccdImageList) {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        // Iteration happens internal to the loop, as we may delete measuredStars from catalog.
        for (MeasuredStarIterator mi = catalog.begin(); mi != catalog.end();) {
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
    for (FittedStarIterator fi = fittedStarList.begin(); fi != fittedStarList.end();) {
        if ((*fi)->getMeasurementCount() == 0) {
            fi = fittedStarList.erase(fi);
        } else {
            ++fi;
        }
    }

    LOGLS_INFO(_log, "Fitted stars after measurement # cut: " << fittedStarList.size());
    LOGLS_INFO(_log, "Total, valid number of Measured stars: " << totalMeasured << ", " << validMeasured);
}

void Associations::normalizeFittedStars() const {
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

    for (auto &fi : fittedStarList) {
        auto measurementCount = fi->getMeasurementCount();
        fi->x /= measurementCount;
        fi->y /= measurementCount;
        fi->getFlux() /= measurementCount;
        fi->getMag() = utils::nanojanskyToABMagnitude(fi->getFlux());
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

void Associations::setFittedStarColors(std::string dicStarListName, std::string color,
                                       double matchCutArcSec) {
    // decode color string in case it is x-y
    size_t pos_minus = color.find('-');
    bool compute_diff = (pos_minus != string::npos);
    std::string c1, c2;
    c1 = color.substr(0, pos_minus);  // if pos_minus == npos, means "up to the end"
    if (compute_diff) c2 = color.substr(pos_minus + 1, string::npos);
    DicStarList cList(dicStarListName);
    if (!cList.HasKey(c1))
        throw(GastroException("Associations::SetFittedstarColors : " + dicStarListName +
                              " misses  a key named \"" + c1 + "\""));
    if (compute_diff && !cList.HasKey(c2))
        throw(GastroException("Associations::SetFittedstarColors : " + dicStarListName +
                              " misses  a key named \"" + c2 + "\""));
    // we associate in some tangent plane. The reference catalog is expressed on the sky,
    // but FittedStar's may be still in this tangent plane.
    BaseStarList &l1 = (BaseStarList &)fittedStarList;
    AstrometryTransformIdentity id;
    TanRaDecToPixel proj(AstrometryTransformLinear(), getCommonTangentPoint());
    // project or not ?
    AstrometryTransform *id_or_proj = &proj;
    if (fittedStarList.inTangentPlaneCoordinates) id_or_proj = &id;
    // The color List is to be projected:
    TStarList projected_cList((BaseStarList &)cList, proj);
    // Associate
    auto starMatchList = listMatchCollect(Fitted2Base(fittedStarList), (const BaseStarList &)projected_cList,
                                          id_or_proj, matchCutArcSec / 3600);

    LOGLS_INFO(_log, "Matched " << starMatchList->size() << '/' << fittedStarList.size()
                                << " FittedStars to color catalog");
    // Evaluate and assign colors.
    for (auto i = starMatchList->begin(); i != starMatchList->end(); ++i) {
        BaseStar *s1 = i->s1;
        FittedStar *fs = dynamic_cast<FittedStar *>(s1);
        BaseStar *s2 = i->s2;
        const TStar *ts = dynamic_cast<const TStar *>(s2);
        const DicStar *ds = dynamic_cast<const DicStar *>(ts->get_original());
        fs->color = ds->getval(c1);
        if (compute_diff) fs->color -= ds->getval(c2);
    }
}

#endif /* TODO */
}  // namespace jointcal
}  // namespace lsst
