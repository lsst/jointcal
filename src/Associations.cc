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

static double sqr(const double &x) {return x * x;}

namespace lsst {
namespace jointcal {

// Source selection is performed in the python, so Associations' constructor is just initializing couple of variables
Associations::Associations()
    : nshoots_(0), nb_photref_associations(0)
{
    commonTangentPoint = Point(0, 0);
}


//bool Associations::AddImage(const std::string &ReducedImageName)
//{
//  ReducedImage ri(ReducedImageName);
//  return AddImage(ri);
//}


bool Associations::AddImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
                            const PTR(lsst::afw::image::TanWcs) wcs,
                            const PTR(lsst::afw::image::VisitInfo) visitInfo,
                            const lsst::afw::geom::Box2I &bbox,
                            const std::string &filter,
                            const PTR(lsst::afw::image::Calib) calib,
                            const int &visit,
                            const int &ccd,
                            const PTR(lsst::jointcal::JointcalControl) control)
{

//  std::cout << " considering image " << ri.Name() << std::endl;

    /* if commonTangentPoint was never initialized
       take the one from this image */
    if (commonTangentPoint.x  == 0. && commonTangentPoint.y == 0.)
    {
        lsst::daf::base::PropertyList::Ptr wcsMeta = wcs->getFitsMetadata();
        double crval1 = wcsMeta->get<double>("CRVAL1");
        double crval2 = wcsMeta->get<double>("CRVAL2");
        commonTangentPoint = Point(crval1, crval2);
        std::cout << "setting commonTangentPoint" << commonTangentPoint << std::endl;
    }

    std::shared_ptr<CcdImage> ccdImage(new CcdImage(Ri, commonTangentPoint, wcs, visitInfo, bbox, filter, calib, visit, ccd, control->sourceFluxField));
//  CcdImage *ccdImage = new CcdImage(Ri, commonTangentPoint, wcs, meta, bbox, filter, calib, visit, ccd, control->sourceFluxField);
    ccdImageList.push_back(ccdImage);
    std::cout << " we have " << ccdImage->WholeCatalog().size()
              << " objects in this catalog " << visit << " " << ccd << std::endl;
    return true;
}

void Associations::AssociateCatalogs(const double MatchCutInArcSec,
                                     const bool UseFittedList,
                                     const bool EnlargeFittedList)
{
    double matchCut = MatchCutInArcSec;

    std::cout << " associating using a cut of " << matchCut << " arcsec" << std::endl;

    // clear catalogs used for previous fits if any
    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        (*i)->CatalogForFit().clear();
    }

    if (!UseFittedList) fittedStarList.clear();
    else // clear measurement counts and associations to refstars.
    {
        for (FittedStarIterator i = fittedStarList.begin();
                i != fittedStarList.end(); ++i)
        {
            (*i)->ClearBeforeAssoc();
        }
    }

    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        CcdImage &ccdImage = **i;

        const Gtransfo *toCommonTangentPlane =
            ccdImage.Pix2CommonTangentPlane();

        /* clear the catalog to fit and copy the whole catalog into it.
        this enables to reassociate from scratch after a fit
            */

        ccdImage.CatalogForFit().clear();
        ccdImage.WholeCatalog().CopyTo(ccdImage.CatalogForFit());
        MeasuredStarList &catalog = ccdImage.CatalogForFit();

        // associate with previous lists
        /* to speed up the match (more precisely the contruction of the
        FastFinder), select in the fittedStarList the objects that
         are within reach of the current ccdImage
            */
        Frame ccdImageFrameCPT =
            ApplyTransfo(ccdImage.ImageFrame(), *toCommonTangentPlane, LargeFrame);
        ccdImageFrameCPT = ccdImageFrameCPT.Rescale(1.10); // add 10 % margin.
        /* we cannot use FittedStarList::ExtractInFrame, because it does an
        actual copy, which we don't want here: we want the pointers in
         the StarMatch to refer to fittedStarList elements. */
        FittedStarList toMatch;

        for (FittedStarCIterator i = fittedStarList.begin();
                i != fittedStarList.end(); ++i)
        {
            if (ccdImageFrameCPT.InFrame(**i))
            {
                toMatch.push_back(*i);
            }
        }


        // divide by 3600 because coordinates in CTP are in degrees.
        StarMatchList *smList = ListMatchCollect(Measured2Base(catalog),
                                Fitted2Base(toMatch),
                                toCommonTangentPlane,
                                matchCut / 3600.);

        /* should check what this RemoveAmbiguities does... */
//      if (Preferences().cleanMatches)
        smList->RemoveAmbiguities(*toCommonTangentPlane);

        /* associate MeasuredStar -> FittedStar using the
        surviving matches */

        int matchedCount = 0;
        for (StarMatchIterator i = smList->begin(); i != smList->end(); ++i)
        {
            StarMatch &starMatch = *i;
            const BaseStar &bs = *starMatch.s1;
            const MeasuredStar &ms_const = dynamic_cast<const MeasuredStar &>(bs);
            MeasuredStar &ms = const_cast<MeasuredStar &>(ms_const);
            const BaseStar &bs2 = *starMatch.s2;
            const FittedStar &fs_const = dynamic_cast<const FittedStar &>(bs2);
            FittedStar &fs = const_cast<FittedStar &>(fs_const);
            ms.SetFittedStar(&fs);
            matchedCount++;

            //      if(  fs->Distance(toCommonTangentPlane->apply(*ms)) > 1.5/3600. )
            //        cout << " OUTLIER: " << ccdImage.Shoot() << "p" << ccdImage.Chip()
            //         << ms->x << "," << ms->y
            //         << endl;

        }
        std::cout << " matched " << matchedCount << " objects"
                  << " in " << ccdImage.Name() << std::endl;
        // delete the matches
        delete smList;

        // add unmatched objets to FittedStarList
        int unMatchedCount = 0;
        for (MeasuredStarIterator i = catalog.begin(); i != catalog.end(); ++i)
        {
            MeasuredStar &mstar = **i;
            // to check if it was matched, just check if it has
            // a fittedStar Pointer assigned
            if (mstar.GetFittedStar()) continue;
            if (EnlargeFittedList)
            {
                FittedStar *fs = new FittedStar(mstar);
                // transform coordinates to CommonTangentPlane
                toCommonTangentPlane->TransformPosAndErrors(*fs, *fs);
                //          fs->Apply(*toCommonTangentPlane);
                fittedStarList.push_back(fs);
                mstar.SetFittedStar(fs);
            }
            unMatchedCount++;
        }
        std::cout << " unmatched objects :" << unMatchedCount << std::endl;
        std::cout << " ************" << std::endl;
    } // end of loop on CcdImage's

    AssignMags();
}
void Associations::CollectRefStars(const bool ProjectOnTP)
{

    // compute the frame on the CTP that contains all input images
    Frame tangentPlaneFrame;

    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        CcdImage &ccdImage = **i;
        Frame CTPFrame = ApplyTransfo(ccdImage.ImageFrame(), *ccdImage.Pix2CommonTangentPlane(), LargeFrame);
        if (tangentPlaneFrame.Area() == 0) tangentPlaneFrame = CTPFrame;
        else tangentPlaneFrame += CTPFrame;
    }

    // DEBUG
    std::cout << " commonTangentPlane Frame " << tangentPlaneFrame
              << " Tangent Point " << commonTangentPoint << std::endl;
    // END DEBUG

    // convert tangent plane coordinates to RaDec:
    GtransfoLin identity;
    TanPix2RaDec CTP2RaDec(identity, commonTangentPoint);
    Frame raDecFrame = ApplyTransfo(tangentPlaneFrame, CTP2RaDec, LargeFrame);

    std::cout << "raDecFrame : " << raDecFrame << std::endl;

    // collect in USNO catalog
    BaseStarList usno;
    UsnoRead(raDecFrame, RColor, usno);
    if (usno.size() == 0)
    {
        throw (LSST_EXCEPT(pex::exceptions::InvalidParameterError, " Nothing read from USNO : stop here "));
    }
    for (BaseStarCIterator i = usno.begin(); i != usno.end(); ++i)
    {
        const BaseStar &b = **i;
        //USNO position uncertainties are set in
        // poloka-core/src/usnoutils.cc (readusno)
        RefStar *r = new RefStar(b, b);
        refStarList.push_back(r);
    }
    // project on CTP (i.e. RaDec2CTP), in degrees
    TanRaDec2Pix RaDec2CTP(identity, commonTangentPoint);

    if (ProjectOnTP)
    {
        refStarList.ApplyTransfo(RaDec2CTP); // also transforms errors in FatPoint
        GtransfoIdentity id;
        AssociateRefStars(usnoMatchCut, &id);
    }
    else
    {
        AssociateRefStars(usnoMatchCut, &RaDec2CTP);
    }
}

void Associations::CollectLSSTRefStars(lsst::afw::table::SortedCatalogT< lsst::afw::table::SimpleRecord > &Ref, std::string filter)
{
    if (Ref.size() == 0)
    {
        throw (LSST_EXCEPT(pex::exceptions::InvalidParameterError, " Reference catalog is empty : stop here "));
    }

//  auto coordKey = Ref.getSchema().find<lsst::afw::coord::Coord>("coord").key;
// Same syntax as the following line but with auto :  auto coordKey = afwTable::CoordKey(Ref.getSchema()["coord"]);
    afw::table::CoordKey coordKey = Ref.getSchema()["coord"];
    auto fluxKey = Ref.getSchema().find<double>(filter + "_flux").key;
//  auto fluxSigmaKey = Ref.getSchema().find<double>(filter + "_fluxSigma").key;

    for (auto i = Ref.begin(); i != Ref.end(); i++)
    {
        lsst::afw::coord::Coord coord = i->get(coordKey);
        double flux = i->get(fluxKey);
//  double fluxErr = i->get(fluxSigmaKey);
        double mag = lsst::afw::image::abMagFromFlux(flux);
        double ra = lsst::afw::geom::radToDeg(coord.getLongitude());
        double dec = lsst::afw::geom::radToDeg(coord.getLatitude());
//  std::cout << flux/fluxErr << " " << mag << std::endl;
//  if (flux/fluxErr < 10.0 || mag > 20. || mag < 16.) {
//      continue;
//  }
        BaseStar s(ra, dec, mag);
        // cook up errors: 100 mas per cooordinate
        s.vx = sqr(0.1 / 3600 / cos(coord.getLatitude()));
        s.vy = sqr(0.1 / 3600);
        s.vxy = 0.;
        RefStar *r = new RefStar(s, s);
        refStarList.push_back(r);
    }

    // project on CTP (i.e. RaDec2CTP), in degrees
    GtransfoLin identity;
    TanRaDec2Pix RaDec2CTP(identity, commonTangentPoint);

    AssociateRefStars(usnoMatchCut, &RaDec2CTP);
}

const lsst::afw::geom::Box2D Associations::GetRaDecBBox()
{
    // compute the frame on the CTP that contains all input images
    Frame tangentPlaneFrame;

    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        CcdImage &ccdImage = **i;
        Frame CTPFrame = ApplyTransfo(ccdImage.ImageFrame(), *ccdImage.Pix2CommonTangentPlane(), LargeFrame);
        if (tangentPlaneFrame.Area() == 0) tangentPlaneFrame = CTPFrame;
        else tangentPlaneFrame += CTPFrame;
    }

    // convert tangent plane coordinates to RaDec:
    GtransfoLin identity;
    TanPix2RaDec CTP2RaDec(identity, commonTangentPoint);
    Frame raDecFrame = ApplyTransfo(tangentPlaneFrame, CTP2RaDec, LargeFrame);

    lsst::afw::geom::Point<double> min(raDecFrame.xMin, raDecFrame.yMin);
    lsst::afw::geom::Point<double> max(raDecFrame.xMax, raDecFrame.yMax);
    lsst::afw::geom::Box2D box(min, max);

    return box;

}

void Associations::AssociateRefStars(const double &MatchCutInArcSec,
                                     const Gtransfo* T)
{
    //DEBUG
    //  std::cout << "****" << std::endl << refStarList << std::endl;
    // std::cout << "fitted list size " << fittedStarList.size() << endl;
    // std::cout << fittedStarList << std::endl;

    // clear previous associations if any :
    for (RefStarIterator i = refStarList.begin(); i != refStarList.end(); ++i)
    {
        (*i)->SetFittedStar(*(FittedStar *) NULL );
    }

    std::cout << " AssociateRefStars : MatchCutInArcSec " << MatchCutInArcSec << std::endl;

    // associate with FittedStars
    // 3600 because coordinates are in degrees (in CTP).

    StarMatchList *smList = ListMatchCollect(Ref2Base(refStarList),
                            Fitted2Base(fittedStarList),
                            T,
                            MatchCutInArcSec / 3600.);

    if (cleanMatches)
    {
        std::cout << " number of matches before removing ambiguities "
                  << smList->size() << std::endl;
        smList->RemoveAmbiguities(*T);
        std::cout << " number of matches after removing ambiguities "
                  << smList->size() << std::endl;
    }

    // actually associate things
    for (StarMatchIterator i = smList->begin(); i != smList->end(); ++i)
    {
        StarMatch &starMatch = *i;
        const BaseStar &bs = *starMatch.s1;
        const RefStar &rs_const = dynamic_cast<const RefStar &>(bs);
        RefStar &rs = const_cast<RefStar &>(rs_const);
        const BaseStar &bs2 = *starMatch.s2;
        const FittedStar &fs_const = dynamic_cast<const FittedStar &>(bs2);
        FittedStar &fs = const_cast<FittedStar &>(fs_const);
        //rs->SetFittedStar(*fs);
        fs.SetRefStar(&rs);
    }

    std::cout << " associated " << smList->size() << " REFERENCE stars "
              << " among a list of " << refStarList.size() << std::endl;
    delete smList;
}
void Associations::SelectFittedStars()
{
    std::cout << " number of possible fitted star before cutting on # of measurements " << fittedStarList.size() << std::endl;
    std::cout << " INFO: min # of measurements " <<  minMeasurementCount << std::endl;
    /* first pass : remove objects that have less than a
       certain number of measurements.
    */
    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        CcdImage &ccdImage = **i;
        MeasuredStarList &catalog = ccdImage.CatalogForFit();
        for (MeasuredStarIterator mi = catalog.begin(); mi != catalog.end(); )
        {
            MeasuredStar &mstar = **mi;

            const FittedStar *fstar = mstar.GetFittedStar();
            if (!fstar) {++mi; continue;}
            int nmes = fstar->MeasurementCount(); // DEBUG

            /*  keep FittedStar's which either have a minimum number of
                measurements, or are matched to a RefStar (i.e. USNO)
            */
            if (!fstar->GetRefStar()
                    &&  fstar->MeasurementCount() < minMeasurementCount)
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

    // just for printouts:
    // fittedStarList.sort(&DecreasingMeasurementCount);

    std::cout
            << " number of possible fitted star after cutting on # of measurements "
            << fittedStarList.size() << std::endl;
}

void Associations::AssignMags()
{
    for (CcdImageIterator i = ccdImageList.begin(); i != ccdImageList.end(); ++i)
    {
        CcdImage &ccdImage = **i;
        MeasuredStarList &catalog = ccdImage.CatalogForFit();
        for (MeasuredStarIterator mi = catalog.begin();
                mi != catalog.end(); ++mi)
        {
            MeasuredStar &mstar = **mi;
            const FittedStar *fstar = mstar.GetFittedStar();
            if (!fstar) continue;
            FittedStar *f = const_cast<FittedStar *>(fstar);
            f->AddMagMeasurement(mstar.Mag(), mstar.MagWeight());
        }
    }
}


void Associations::DeprojectFittedStars()
{
    /* by default, Associations::fittedStarList is expressed on the
       Associations::commonTangentPlane. For AstromFit, we need it on
       the sky */
    if (!fittedStarList.inTangentPlaneCoordinates)
    {
        std::cout << "WARNING: Associations::DeprojectFittedStars : Fitted stars are already in sidereal coordinates, nothing done " << std::endl;
        return;
    }

    TanPix2RaDec ctp2Sky(GtransfoLin(), CommonTangentPoint());
    fittedStarList.ApplyTransfo(ctp2Sky);
    fittedStarList.inTangentPlaneCoordinates = false;
}


#ifdef STORAGE
void Associations::CollectMCStars(int realization)
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
        MeasuredStarList& catalog = ccdImage.CatalogForFit();

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


//void Associations::CheckMCStars()
//{
//  CcdImageIterator I;
//  StarMatchIterator smI;
//
//  for(I=ccdImageList.begin();I!=ccdImageList.end();I++)
//    {
//      CcdImage& ccdImage = **I;
//      string dbimdir = ccdImage.Dir();
//      string mcstars = dbimdir + "/mcstars.list";
//      string mctruth = dbimdir + "/mctruth.list";
//      GtransfoIdentity gti;
//      MeasuredStarList& catalog = ccdImage.CatalogForFit();
//
//      MeasuredStarIterator msII;
//      for(msII=catalog.begin();msII!=catalog.end();msII++)
//          {
//            MeasuredStar& ms = **msII;
//            cout << " *** >>> "
//                 << ms.GetMCMeas() << " "
//                 << ms.GetMCTruth() << endl;
//          }
//    }
//}


static double getflux(DicStar const& ds, string const& tag)
{
    double rflux = ds.getval(tag);
    if (rflux <= 0.) return -1.;
    return 1.E10 * pow(10., -0.4 * rflux);
}

void Associations::CollectPhotometricRefStars(string const& catalogname)
{
    cout << "[Associations::CollectPhotometricRefStars()] ... ";
    photRefStarList.clear();

    DicStarList dsl(catalogname);
    DicStarIterator I;

    double rflux;
    vector<double> refflux;
    refflux.resize(10);

    for (I = dsl.begin(); I != dsl.end(); I++)
    {
        DicStar& ds = **I;

        // I am here (debugging the star positions)
        //       cout << "[DEBUG] x,y="
        //       << ds.x << " " << ds.y
        //       << endl;

        refflux[0] = getflux(ds, "mu");
        refflux[1] = getflux(ds, "mg");
        refflux[2] = getflux(ds, "mr");
        refflux[3] = getflux(ds, "mi");
        refflux[4] = getflux(ds, "mz");

        refflux[5] = getflux(ds, "U");
        refflux[6] = getflux(ds, "B");
        refflux[7] = getflux(ds, "V");
        refflux[8] = getflux(ds, "R");
        refflux[9] = getflux(ds, "I");

        RefStar* rs = new RefStar(ds, ds);
        rs->AssignRefFluxes(refflux);

        unsigned int index = (unsigned int)ds.getval("index");
        rs->Index() = index;

        photRefStarList.push_back(rs);
    }

    cout << photRefStarList.size() << " stars read" << endl;

    // Now, If I understand well, all the FittedStar coordinates
    // are in CTP coordinates ==> we have to transform our
    // (RA,DEC) coordinates into CTP coords...
    GtransfoLin identity;
    TanRaDec2Pix RaDec2CTP(identity, commonTangentPoint);
    photRefStarList.ApplyTransfo(RaDec2CTP);

    cout << "[DEBUG DEBUB DEBUG]" << endl;
    RaDec2CTP.dump(cout);
    cout << "[DEBUG DEBUB DEBUG]" << endl;
}


void Associations::AssociatePhotometricRefStars(double MatchCutInArcSec)
{
    cout << "[Associations::AssociatePhotometricRefStars]" << endl;
    cout << "[DEBUG]  -> nb fs=" << fittedStarList.size()
         << " nb refstars="  << photRefStarList.size() << endl;

    // clear previous associations
    RefStarIterator I;
    for (I = photRefStarList.begin(); I != photRefStarList.end(); I++)
    {
        (*I)->SetFittedStar( *(FittedStar*)NULL ); // mais c'est horrible (!)
        //      cout << " XTP,YTP=" << (*I)->x << " " << (*I)->y << endl;
    }

    GtransfoIdentity gtransfoIdentity;
    StarMatchList* smList = ListMatchCollect(Ref2Base(photRefStarList),
                            Fitted2Base(fittedStarList),
                            &gtransfoIdentity,
                            MatchCutInArcSec / 3600.);

    if (Preferences().cleanMatches)
    {
        cout << "  (*) removing ambiguities: "
             << smList->size();
        smList->RemoveAmbiguities(gtransfoIdentity);
        cout << " --> "  << smList->size() << endl;
    }

    // associate the RefStars and the fitted stars ...
    StarMatchIterator smI;
    for (smI = smList->begin(); smI != smList->end(); smI++)
    {
        StarMatch& sm = *smI; // et pas **smI cette fois-ci. C'est d'un penible !
        BaseStar* bs = sm.s1;
        RefStar* rs = dynamic_cast<RefStar*>(bs); // et dynamic_cast<RefStar*>(sm.s1) ne compile pas
        bs = sm.s2;
        FittedStar* fs = dynamic_cast<FittedStar*>(bs);
        fs->SetRefStar(*rs);
    }

    cout << "  (*) done. Associated " << smList->size() << " STANDARD stars "
         << "among a list of " << photRefStarList.size() << " objects."
         << endl;

    nb_photref_associations = smList->size();

    delete smList;
}



//void Associations::SetRefPhotFactor(int chip, double photfact)
//{
//
//}



#include "tstar.h"

void Associations::SetFittedStarColors(std::string DicStarListName,
                                       std::string Color,
                                       const double &MatchCutArcSec)
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
    TanRaDec2Pix proj(GtransfoLin(), CommonTangentPoint());
    // project or not ?
    Gtransfo *id_or_proj = &proj;
    if (fittedStarList.inTangentPlaneCoordinates) id_or_proj = &id;
    // The color List is to be projected:
    TStarList projected_cList((BaseStarList &) cList, proj);
    // Associate
    StarMatchList *sm = ListMatchCollect(Fitted2Base(fittedStarList),
                                         (const BaseStarList &) projected_cList,
                                         id_or_proj,
                                         MatchCutArcSec / 3600);

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

void Associations::CollectMCStars(int realization)
{
    CcdImageIterator I;
    StarMatchIterator smI;

    for (I = ccdImageList.begin(); I != ccdImageList.end(); I++)
    {
        CcdImage& ccdImage = **I;
        string dbimdir = ccdImage.Dir();
        string mcstars = dbimdir + "/mc/mcstars.list";
        string mctruth = dbimdir + "/mc/mctruth.list";

        if (realization >= 0)
        {
            stringstream sstrm;
            sstrm << dbimdir << "/mc/mcstars_" << realization << ".list";
            mcstars = sstrm.str();
            sstrm.str("");
            sstrm << dbimdir << "/mc/mctruth_" << realization << ".list";
            mctruth = sstrm.str();
        }

        GtransfoIdentity gti;
        MeasuredStarList& catalog = ccdImage.CatalogForFit();

        BaseStarWithErrorList mcstarlist(mcstars);
        StarMatchList *smList = ListMatchCollect(Measured2Base(catalog),
                                BaseStarWithError2Base(mcstarlist),
                                &gti, 1. /*pixel*/);
        if (smList)
            for (smI = smList->begin(); smI != smList->end(); smI++)
            {
                StarMatch& sm = *smI;
                BaseStar* bs = sm.s1;
                MeasuredStar* mstar = dynamic_cast<MeasuredStar* >(bs);
                bs = sm.s2;
                BaseStarWithError* mcstar = dynamic_cast<BaseStarWithError* >(bs);
                mstar->SetMCMeas(mcstar);
            }
        else
            cout << "[Associations::CollectMCStars] Unable to match MCStars w/ catalog !" << endl;
        delete smList;

        BaseStarWithErrorList mctruthlist(mctruth);
        //      DictStarList mctruthlist(mctruth);
        smList = ListMatchCollect(Measured2Base(catalog),
                                  BaseStarWithError2Base(mctruthlist),
                                  &gti, 1. /* pixel ? */);
        if (smList)
            for (smI = smList->begin(); smI != smList->end(); smI++)
            {
                StarMatch& sm = *smI;
                BaseStar* bs = sm.s1;
                MeasuredStar* mstar = dynamic_cast<MeasuredStar*>(bs);
                bs = sm.s2;
                BaseStarWithError* mcstar = dynamic_cast<BaseStarWithError* >(bs);
                mstar->SetMCTruth(mcstar);
            }
        else
            cout << "[Associations::CollectMCStars] Unable to match MCTruth w/ catalog !" << endl;
        delete smList;

        //      MeasuredStarIterator msII;
        //      for(msII=catalog.begin();msII!=catalog.end();msII++)
        //    {
        //      MeasuredStar& ms = **msII;
        //      cout << " *** >>> "
        //           << ms.GetMCMeas() << " "
        //           << ms.GetMCTruth() << endl;
        //    }

    }
}
#endif /* STORAGE */

}
} // end of namespaces
