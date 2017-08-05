// -*- C++ -*-
#include <iostream>
#include <limits>
#include <sstream>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/ListMatch.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/image/Calib.h"

namespace jointcal = lsst::jointcal;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.Associations");
}

// TODO: Remove this once RFC-356 is implemented and all refcats give fluxes in Maggies.
const double JanskyToMaggy = 3631.0;

namespace lsst {
namespace jointcal {

void Associations::addImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &catalog,
                            std::shared_ptr<lsst::afw::image::TanWcs> wcs,
                            std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo,
                            lsst::afw::geom::Box2I const &bbox, std::string const &filter,
                            std::shared_ptr<afw::image::PhotoCalib> photoCalib,
                            std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
                            std::shared_ptr<lsst::jointcal::JointcalControl> control) {
    auto ccdImage = std::make_shared<CcdImage>(catalog, wcs, visitInfo, bbox, filter, photoCalib, detector,
                                               visit, ccd, control->sourceFluxField);
    ccdImageList.push_back(ccdImage);
    LOGLS_DEBUG(_log, "Catalog " << ccdImage->getName() << " has " << ccdImage->getWholeCatalog().size()
                                 << " objects.");
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
        const Gtransfo *toCommonTangentPlane = ccdImage->getPix2CommonTangentPlane();

        /* clear the catalog to fit and copy the whole catalog into it.
        this allows reassociating from scratch after a fit. */

        ccdImage->getCatalogForFit().clear();
        ccdImage->getWholeCatalog().copyTo(ccdImage->getCatalogForFit());
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();

        // associate with previous lists
        /* to speed up the match (more precisely the contruction of the
        FastFinder), select in the fittedStarList the objects that
         are within reach of the current ccdImage
            */
        Frame ccdImageFrameCPT = applyTransfo(ccdImage->getImageFrame(), *toCommonTangentPlane, LargeFrame);
        ccdImageFrameCPT = ccdImageFrameCPT.rescale(1.10);  // add 10 % margin.
        /* we cannot use FittedStarList::ExtractInFrame, because it does an
        actual copy, which we don't want here: we want the pointers in
         the StarMatch to refer to fittedStarList elements. */
        FittedStarList toMatch;

        for (auto const &fittedStar : fittedStarList) {
            if (ccdImageFrameCPT.inFrame(*fittedStar)) {
                toMatch.push_back(fittedStar);
            }
        }

        // divide by 3600 because coordinates in CTP are in degrees.
        auto starMatchList = listMatchCollect(Measured2Base(catalog), Fitted2Base(toMatch),
                                              toCommonTangentPlane, matchCutInArcSec / 3600.);

        /* should check what this removeAmbiguities does... */
        LOGLS_DEBUG(_log, "Measured-to-Fitted matches before removing ambiguities " << starMatchList->size());
        starMatchList->removeAmbiguities(*toCommonTangentPlane);
        LOGLS_DEBUG(_log, "Measured-to-Fitted matches after removing ambiguities " << starMatchList->size());

        /* associate MeasuredStar -> FittedStar using the
        surviving matches */

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

    assignMags();
}

void Associations::collectRefStars(lsst::afw::table::SortedCatalogT<lsst::afw::table::SimpleRecord> &refCat,
                                   afw::geom::Angle matchCut, std::string const &fluxField,
                                   std::map<std::string, std::vector<double>> const &refFluxMap,
                                   std::map<std::string, std::vector<double>> const &refFluxErrMap) {
    if (refCat.size() == 0) {
        throw(LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          " reference catalog is empty : stop here "));
    }

    afw::table::CoordKey coordKey = refCat.getSchema()["coord"];
    auto fluxKey = refCat.getSchema().find<double>(fluxField).key;
    // Don't blow up if the reference catalog doesn't contain errors.
    afw::table::Key<double> fluxErrKey;
    try {
        fluxErrKey = refCat.getSchema().find<double>(fluxField + "Sigma").key;
    } catch (pex::exceptions::NotFoundError &) {
        LOGLS_WARN(_log, "Flux error field ("
                                 << fluxField << "Sigma"
                                 << ") not found in reference catalog. Not using ref flux errors.");
    }
    std::cout << "Error key: " << fluxErrKey << " valid: " << fluxErrKey.isValid() << std::endl;
    _filterMap.clear();
    _filterMap.reserve(refFluxMap.size());
    size_t nFilters = 0;
    for (auto const &filter : refFluxMap) {
        _filterMap[filter.first] = nFilters;
        nFilters++;
    }

    refStarList.clear();
    for (size_t i = 0; i < refCat.size(); i++) {
        auto const &record = refCat.get(i);

        afw::coord::Coord coord = record->get(coordKey);
        double defaultFlux = record->get(fluxKey) / JanskyToMaggy;
        double defaultFluxErr;
        if (fluxErrKey.isValid()) {
            defaultFluxErr = record->get(fluxErrKey) / JanskyToMaggy;
        } else {
            defaultFluxErr = std::numeric_limits<double>::quiet_NaN();
        }
        std::vector<double> fluxList(nFilters);
        std::vector<double> fluxErrList(nFilters);
        for (auto const &filter : _filterMap) {
            fluxList[filter.second] = refFluxMap.at(filter.first).at(i) / JanskyToMaggy;
            fluxErrList[filter.second] = refFluxErrMap.at(filter.first).at(i) / JanskyToMaggy;
        }
        double ra = lsst::afw::geom::radToDeg(coord.getLongitude());
        double dec = lsst::afw::geom::radToDeg(coord.getLatitude());
        auto star = std::make_shared<RefStar>(ra, dec, defaultFlux, defaultFluxErr, fluxList, fluxErrList);

        // TODO DM-10826: RefCats aren't guaranteed to have position errors.
        // TODO: Need to devise a way to check whether the refCat has position errors
        // TODO: and use them instead, if available.
        // cook up errors: 100 mas per cooordinate
        star->vx = std::pow(0.1 / 3600 / cos(coord.getLatitude()), 2);
        star->vy = std::pow(0.1 / 3600, 2);
        star->vxy = 0.;

        refStarList.push_back(star);
    }

    // project on CTP (i.e. RaDec2CTP), in degrees
    GtransfoLin identity;
    TanRaDec2Pix raDec2CTP(identity, _commonTangentPoint);

    associateRefStars(matchCut.asArcseconds(), &raDec2CTP);
}

const lsst::afw::geom::Box2D Associations::getRaDecBBox() {
    // compute the frame on the CTP that contains all input images
    Frame tangentPlaneFrame;

    for (auto const &ccdImage : ccdImageList) {
        Frame CTPFrame =
                applyTransfo(ccdImage->getImageFrame(), *(ccdImage->getPix2CommonTangentPlane()), LargeFrame);
        if (tangentPlaneFrame.getArea() == 0)
            tangentPlaneFrame = CTPFrame;
        else
            tangentPlaneFrame += CTPFrame;
    }

    // convert tangent plane coordinates to RaDec:
    GtransfoLin identity;
    TanPix2RaDec CTP2RaDec(identity, _commonTangentPoint);
    Frame raDecFrame = applyTransfo(tangentPlaneFrame, CTP2RaDec, LargeFrame);

    lsst::afw::geom::Point<double> min(raDecFrame.xMin, raDecFrame.yMin);
    lsst::afw::geom::Point<double> max(raDecFrame.xMax, raDecFrame.yMax);
    lsst::afw::geom::Box2D box(min, max);

    return box;
}

void Associations::associateRefStars(double matchCutInArcSec, const Gtransfo *gtransfo) {
    // associate with FittedStars
    // 3600 because coordinates are in degrees (in CTP).
    auto starMatchList = listMatchCollect(Ref2Base(refStarList), Fitted2Base(fittedStarList), gtransfo,
                                          matchCutInArcSec / 3600.);

    LOGLS_DEBUG(_log, "Refcat matches before removing ambiguities " << starMatchList->size());
    starMatchList->removeAmbiguities(*gtransfo);
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

void Associations::selectFittedStars(int minMeasurements) {
    LOGLS_INFO(_log, "Fitted stars before measurement # cut: " << fittedStarList.size());
    /* first pass : remove objects that have less than a
       certain number of measurements.
    */
    for (auto const &ccdImage : ccdImageList) {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (MeasuredStarIterator mi = catalog.begin(); mi != catalog.end();) {
            MeasuredStar &mstar = **mi;

            auto fstar = mstar.getFittedStar();
            if (!fstar) {
                ++mi;
                continue;
            }

            /*  keep FittedStar's which either have a minimum number of
                measurements, or are matched to a RefStar
            */
            if (!fstar->getRefStar() && fstar->getMeasurementCount() < minMeasurements) {
                auto f = std::const_pointer_cast<FittedStar>(fstar);
                f->getMeasurementCount()--;
                mi = catalog.erase(mi);
            } else
                ++mi;
        }  // end loop on objects in catalog
    }      // end loop on catalogs

    /* now FittedStars with less than minMeasurements should have
       zero measurementCount; */

    for (FittedStarIterator fi = fittedStarList.begin(); fi != fittedStarList.end();) {
        if ((*fi)->getMeasurementCount() == 0)
            fi = fittedStarList.erase(fi);
        else
            ++fi;
    }

    LOGLS_INFO(_log, "Fitted stars after measurement # cut: " << fittedStarList.size());
}

void Associations::assignMags() {
    for (auto const &ccdImage : ccdImageList) {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (auto const &mstar : catalog) {
            auto fstar = mstar->getFittedStar();
            if (!fstar) continue;
            auto f = std::const_pointer_cast<FittedStar>(fstar);
            f->addMagMeasurement(mstar->getMag(), mstar->getMagWeight());
        }
    }
}

void Associations::deprojectFittedStars() {
    /* by default, Associations::fittedStarList is expressed on the
       Associations::commonTangentPlane. For AstrometryFit, we need it on
       the sky */
    if (!fittedStarList.inTangentPlaneCoordinates) {
        LOGLS_WARN(_log,
                   "DeprojectFittedStars: Fitted stars are already in sidereal coordinates, nothing done ");
        return;
    }

    TanPix2RaDec ctp2Sky(GtransfoLin(), getCommonTangentPoint());
    fittedStarList.applyTransfo(ctp2Sky);
    fittedStarList.inTangentPlaneCoordinates = false;
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

        GtransfoIdentity gti;
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
    GtransfoIdentity id;
    TanRaDec2Pix proj(GtransfoLin(), getCommonTangentPoint());
    // project or not ?
    Gtransfo *id_or_proj = &proj;
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

#endif /* STORAGE */
}  // namespace jointcal
}  // namespace lsst
