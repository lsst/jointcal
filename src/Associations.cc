// -*- C++ -*-
//
#include <iostream>
#include <sstream>

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

// TODO: propagate those into python:
const double usnoMatchCut = 3;
const bool cleanMatches = true;
const int minMeasurementCount = 2;

namespace jointcal = lsst::jointcal;

static double sqr(double x) {return x*x;}

namespace lsst {
namespace jointcal {

Associations::Associations()
{
    _commonTangentPoint = Point(0, 0);
}

void Associations::addImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &catalog,
                            std::shared_ptr<lsst::afw::image::TanWcs> wcs,
                            std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo,
                            lsst::afw::geom::Box2I const &bbox,
                            std::string const &filter,
                            std::shared_ptr<lsst::afw::image::Calib> calib,
                            int visit,
                            int ccd,
                            std::shared_ptr<lsst::jointcal::JointcalControl> control)
{
    std::shared_ptr<CcdImage> ccdImage(new CcdImage(catalog, wcs, visitInfo, bbox, filter, calib, visit, ccd, control->sourceFluxField));
    ccdImageList.push_back(ccdImage);
    std::cout << " we have " << ccdImage->getWholeCatalog().size()
              << " objects in this catalog " << visit << " " << ccd << std::endl;
}

void Associations::setCommonTangentPoint(lsst::afw::geom::Point2D const &commonTangentPoint)
{
    _commonTangentPoint = Point(commonTangentPoint.getX(), commonTangentPoint.getY()); // a jointcal::Point
    for (auto &ccdImage: ccdImageList)
        ccdImage->setCommonTangentPoint(_commonTangentPoint);
}

void Associations::associateCatalogs(const double matchCutInArcSec,
                                     const bool UseFittedList,
                                     const bool EnlargeFittedList)
{
    std::cout << " associating using a cut of " << matchCutInArcSec << " arcsec" << std::endl;

    // clear reference stars
    for (auto &item: refStarList)
    {
        item->setFittedStar(nullptr);
    }
    refStarList.clear();

    // clear measurement counts and associations to refstars, but keep fittedStars themselves.
    for (auto &item: fittedStarList)
    {
        item->clearBeforeAssoc();
    }
    // clear fitted stars
    if (!UseFittedList)
        fittedStarList.clear();

    for (auto &ccdImage: ccdImageList)
    {
        const Gtransfo *toCommonTangentPlane = ccdImage->Pix2CommonTangentPlane();

        /* clear the catalog to fit and copy the whole catalog into it.
        this allows reassociating from scratch after a fit. */

        ccdImage->getCatalogForFit().clear();
        ccdImage->getWholeCatalog().CopyTo(ccdImage->getCatalogForFit());
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();

        // associate with previous lists
        /* to speed up the match (more precisely the contruction of the
        FastFinder), select in the fittedStarList the objects that
         are within reach of the current ccdImage
            */
        Frame ccdImageFrameCPT =
            ApplyTransfo(ccdImage->ImageFrame(), *toCommonTangentPlane, LargeFrame);
        ccdImageFrameCPT = ccdImageFrameCPT.Rescale(1.10); // add 10 % margin.
        /* we cannot use FittedStarList::ExtractInFrame, because it does an
        actual copy, which we don't want here: we want the pointers in
         the StarMatch to refer to fittedStarList elements. */
        FittedStarList toMatch;

        for (auto const &fittedStar: fittedStarList)
        {
            if (ccdImageFrameCPT.InFrame(*fittedStar))
            {
                toMatch.push_back(fittedStar);
            }
        }

        // divide by 3600 because coordinates in CTP are in degrees.
        StarMatchList *smList = ListMatchCollect(Measured2Base(catalog),
                                                 Fitted2Base(toMatch),
                                                 toCommonTangentPlane,
                                                 matchCutInArcSec/3600.);

        /* should check what this RemoveAmbiguities does... */
//      if (Preferences().cleanMatches)
        smList->RemoveAmbiguities(*toCommonTangentPlane);

        /* associate MeasuredStar -> FittedStar using the
        surviving matches */

        int matchedCount = 0;
        for (auto const &starMatch: *smList)
        {
            const BaseStar &bs = *starMatch.s1;
            const MeasuredStar &ms_const = dynamic_cast<const MeasuredStar &>(bs);
            MeasuredStar &ms = const_cast<MeasuredStar &>(ms_const);
            const BaseStar &bs2 = *starMatch.s2;
            const FittedStar &fs_const = dynamic_cast<const FittedStar &>(bs2);
            FittedStar &fs = const_cast<FittedStar &>(fs_const);
            ms.SetFittedStar(&fs);
            matchedCount++;
        }
        std::cout << " matched " << matchedCount << " objects"
                  << " in " << ccdImage->getName() << std::endl;
        // delete the matches
        delete smList;

        // add unmatched objets to FittedStarList
        int unMatchedCount = 0;
        for (auto const &mstar: catalog)
        {
            // to check if it was matched, just check if it has
            // a fittedStar Pointer assigned
            if (mstar->GetFittedStar()) continue;
            if (EnlargeFittedList)
            {
                FittedStar *fs = new FittedStar(*mstar);
                // transform coordinates to CommonTangentPlane
                toCommonTangentPlane->TransformPosAndErrors(*fs, *fs);
                fittedStarList.push_back(fs);
                mstar->SetFittedStar(fs);
            }
            unMatchedCount++;
        }
        std::cout << " unmatched objects :" << unMatchedCount << std::endl;
        std::cout << " ************" << std::endl;
    } // end of loop on CcdImage's

    assignMags();
}

void Associations::collectRefStars(lsst::afw::table::SortedCatalogT< lsst::afw::table::SimpleRecord > &refCat,
                                   std::string const &fluxField)
{
    if (refCat.size() == 0)
    {
        throw (LSST_EXCEPT(pex::exceptions::InvalidParameterError, " reference catalog is empty : stop here "));
    }

//  auto coordKey = refCat.getSchema().find<lsst::afw::coord::Coord>("coord").key;
// Same syntax as the following line but with auto :  auto coordKey = afwTable::CoordKey(refCat.getSchema()["coord"]);
    afw::table::CoordKey coordKey = refCat.getSchema()["coord"];
    auto fluxKey = refCat.getSchema().find<double>(fluxField).key;

    for (auto const &i: refCat)
    {
        lsst::afw::coord::Coord coord = i.get(coordKey);
        double flux = i.get(fluxKey);
        double mag = lsst::afw::image::abMagFromFlux(flux);
        double ra = lsst::afw::geom::radToDeg(coord.getLongitude());
        double dec = lsst::afw::geom::radToDeg(coord.getLatitude());
        BaseStar s(ra, dec, mag);
        // cook up errors: 100 mas per cooordinate

        // TODO: What is this? Why are we making fake errors here?

        s.vx = sqr(0.1/3600/cos(coord.getLatitude()));
        s.vy = sqr(0.1/3600);
        s.vxy = 0.;
        RefStar *r = new RefStar(s);
        refStarList.push_back(r);
    }

    // project on CTP (i.e. RaDec2CTP), in degrees
    GtransfoLin identity;
    TanRaDec2Pix raDec2CTP(identity, _commonTangentPoint);

    associateRefStars(usnoMatchCut, &raDec2CTP);
}

const lsst::afw::geom::Box2D Associations::getRaDecBBox()
{
    // compute the frame on the CTP that contains all input images
    Frame tangentPlaneFrame;

    for (auto const &ccdImage: ccdImageList)
    {
        Frame CTPFrame = ApplyTransfo(ccdImage->ImageFrame(), *(ccdImage->Pix2CommonTangentPlane()), LargeFrame);
        if (tangentPlaneFrame.Area() == 0) tangentPlaneFrame = CTPFrame;
        else tangentPlaneFrame += CTPFrame;
    }

    // convert tangent plane coordinates to RaDec:
    GtransfoLin identity;
    TanPix2RaDec CTP2RaDec(identity, _commonTangentPoint);
    Frame raDecFrame = ApplyTransfo(tangentPlaneFrame, CTP2RaDec, LargeFrame);

    lsst::afw::geom::Point<double> min(raDecFrame.xMin, raDecFrame.yMin);
    lsst::afw::geom::Point<double> max(raDecFrame.xMax, raDecFrame.yMax);
    lsst::afw::geom::Box2D box(min, max);

    return box;

}

void Associations::associateRefStars(double matchCutInArcSec, const Gtransfo* gtransfo)
{
    // associate with FittedStars
    // 3600 because coordinates are in degrees (in CTP).
    StarMatchList *smList = ListMatchCollect(Ref2Base(refStarList),
                                             Fitted2Base(fittedStarList),
                                             gtransfo,
                                             matchCutInArcSec/3600.);

    if (cleanMatches)
    {
        std::cout << " number of refcat matches before removing ambiguities "
                  << smList->size() << std::endl;
        smList->RemoveAmbiguities(*gtransfo);
        std::cout << " number of refcat matches after removing ambiguities "
                  << smList->size() << std::endl;
    }

    // actually associate things
    for (auto const &starMatch: *smList)
    {
        const BaseStar &bs = *starMatch.s1;
        const RefStar &rs_const = dynamic_cast<const RefStar &>(bs);
        RefStar &rs = const_cast<RefStar &>(rs_const);
        const BaseStar &bs2 = *starMatch.s2;
        const FittedStar &fs_const = dynamic_cast<const FittedStar &>(bs2);
        FittedStar &fs = const_cast<FittedStar &>(fs_const);
        //rs->SetFittedStar(*fs);
        fs.setRefStar(&rs);
    }

    std::cout << " associated " << smList->size() << " reference stars "
              << " among a list of " << refStarList.size() << std::endl;
    delete smList;
}

void Associations::selectFittedStars()
{
    std::cout << " number of possible fitted star before cutting on # of measurements " << fittedStarList.size() << std::endl;
    std::cout << " INFO: min # of measurements " <<  minMeasurementCount << std::endl;
    /* first pass : remove objects that have less than a
       certain number of measurements.
    */
    for (auto const &ccdImage: ccdImageList)
    {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (MeasuredStarIterator mi = catalog.begin(); mi != catalog.end(); )
        {
            MeasuredStar &mstar = **mi;

            const FittedStar *fstar = mstar.GetFittedStar();
            if (!fstar) {++mi; continue;}
            int nmes = fstar->MeasurementCount(); // DEBUG

            /*  keep FittedStar's which either have a minimum number of
                measurements, or are matched to a RefStar
            */
            if (!fstar->getRefStar() &&  fstar->MeasurementCount() < minMeasurementCount)
            {
                FittedStar *f = const_cast<FittedStar *>(fstar);
                f->MeasurementCount()--;
                mi = catalog.erase(mi);
                // DEBUG
                if (fstar && fstar->MeasurementCount() != nmes - 1)
                {
                    std::cout << " ca craint " << std::endl;
                }
            }
            else ++mi;
        }// end loop on objects in catalog
    } // end loop on catalogs

    /* now FittedStars with less than minMeasurementCount should have
       zero MeasurementCount(); */

    for (FittedStarIterator fi = fittedStarList.begin();
            fi != fittedStarList.end();  )
    {
        if ((*fi)->MeasurementCount() == 0) fi = fittedStarList.erase(fi);
        else ++fi;
    }

    std::cout
            << " number of possible fitted star after cutting on # of measurements "
            << fittedStarList.size() << std::endl;
}

void Associations::assignMags()
{
    for (auto const &ccdImage: ccdImageList)
    {
        MeasuredStarList &catalog = ccdImage->getCatalogForFit();
        for (auto const &mstar: catalog)
        {
            const FittedStar *fstar = mstar->GetFittedStar();
            if (!fstar) continue;
            FittedStar *f = const_cast<FittedStar *>(fstar);
            f->AddMagMeasurement(mstar->Mag(), mstar->MagWeight());
        }
    }
}

void Associations::deprojectFittedStars()
{
    /* by default, Associations::fittedStarList is expressed on the
       Associations::commonTangentPlane. For AstromFit, we need it on
       the sky */
    if (!fittedStarList.inTangentPlaneCoordinates)
    {
        std::cout << "WARNING: Associations::DeprojectFittedStars : Fitted stars are already in sidereal coordinates, nothing done " << std::endl;
        return;
    }

    TanPix2RaDec ctp2Sky(GtransfoLin(), getCommonTangentPoint());
    fittedStarList.ApplyTransfo(ctp2Sky);
    fittedStarList.inTangentPlaneCoordinates = false;
}

#ifdef STORAGE
void Associations::collectMCStars(int realization)
{
    cout << "[Associations::CollectMCStars]" << endl;
    CcdImageIterator I;
    StarMatchIterator smI;

    for (I = ccdImageList.begin(); I != ccdImageList.end(); I++)
    {
        CcdImage& ccdImage = **I;
        string dbimdir = ccdImage.Dir();
        string mctruth = dbimdir + "/mc/mctruth.list";

        if (realization >= 0)
        {
            stringstream sstrm;
            sstrm << dbimdir << "/mc/mctruth_" << realization << ".list";
            mctruth = sstrm.str();
        }

        GtransfoIdentity gti;
        MeasuredStarList& catalog = ccdImage.getCatalogForFit();

        //      BaseStarWithErrorList mctruthlist(mctruth);
        DicStarList mctruthlist(mctruth);
        StarMatchList* smList = ListMatchCollect(Measured2Base(catalog),
                                Dic2Base(mctruthlist),
                                &gti, 1. /* pixel ? */);
        if (smList)
            for (smI = smList->begin(); smI != smList->end(); smI++)
            {
                StarMatch& sm = *smI;
                BaseStar* bs = sm.s1;
                MeasuredStar* mstar = dynamic_cast<MeasuredStar*>(bs);
                bs = sm.s2;
                DicStar* dstar = dynamic_cast<DicStar* >(bs);
                BaseStarWithError* mcstar = new BaseStarWithError(*bs);
                mcstar->GetMCInfo().iflux = dstar->getval("iflux");
                mcstar->GetMCInfo().tflux = dstar->getval("sflux");
                /*
                mstar->SetMCTruth(mcstar);
                mstar->SetMCMeas(mcstar);
                */
            }
        else
            cout << "[Associations::CollectMCStars] Unable to match MCTruth w/ catalog !" << endl;
        delete smList;
    }
}

static double getflux(DicStar const& ds, std::string const& tag)
{
    double rflux = ds.getval(tag);
    if (rflux <= 0.) return -1.;
    return 1.E10 * pow(10., -0.4*rflux);
}

void Associations::setFittedStarColors(std::string DicStarListName,
                                       std::string Color,
                                       double MatchCutArcSec)
{
    // decode color string in case it is x-y
    size_t pos_minus = Color.find('-');
    bool compute_diff = (pos_minus != string::npos);
    std::string c1, c2;
    c1 = Color.substr(0, pos_minus); // if pos_minus == npos, means "up to the end"
    if (compute_diff)
        c2 = Color.substr(pos_minus + 1, string::npos);
    DicStarList cList(DicStarListName);
    if (!cList.HasKey(c1)) throw(GastroException("Associations::SetFittedstarColors : " + DicStarListName + " misses  a key named \"" + c1 + "\""));
    if (compute_diff && !cList.HasKey(c2))
        throw(GastroException("Associations::SetFittedstarColors : " + DicStarListName + " misses  a key named \"" + c2 + "\""));
    // we associate in some tangent plane. The reference catalog is expressed on the sky,
    // but FittedStar's may be still in this tangent plane.
    BaseStarList &l1 = (BaseStarList &) fittedStarList;
    GtransfoIdentity id;
    TanRaDec2Pix proj(GtransfoLin(), getCommonTangentPoint());
    // project or not ?
    Gtransfo *id_or_proj = &proj;
    if (fittedStarList.inTangentPlaneCoordinates) id_or_proj = &id;
    // The color List is to be projected:
    TStarList projected_cList((BaseStarList &) cList, proj);
    // Associate
    StarMatchList *sm = ListMatchCollect(Fitted2Base(fittedStarList),
                                         (const BaseStarList &) projected_cList,
                                         id_or_proj,
                                         MatchCutArcSec/3600);

    cout << "INFO : matched " << sm->size() << '/' << fittedStarList.size()
         << " FittedStars to color catalog" << endl;
    // Evaluate and assign colors.
    for (auto i = sm->begin(); i != sm->end(); ++i)
    {
        BaseStar *s1 = i->s1;
        FittedStar *fs = dynamic_cast<FittedStar*>(s1);
        BaseStar *s2 = i->s2;
        const TStar* ts = dynamic_cast<const TStar*>(s2);
        const DicStar *ds = dynamic_cast<const DicStar*>(ts->get_original());
        fs->color = ds->getval(c1);
        if (compute_diff) fs->color -= ds->getval(c2);
    }
}

#endif /* STORAGE */

}
} // end of namespaces
