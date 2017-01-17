#include <assert.h>
#include <string>
#include <sstream>
#include <math.h>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/afw/image/Image.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Point.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Angle.h"

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;

namespace lsst {
namespace jointcal {

static double sq(double x) { return x*x;}

void CcdImage::LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat, const std::string &fluxField)
{
    auto xKey = Cat.getSchema().find<double>("slot_Centroid_x").key;
    auto yKey = Cat.getSchema().find<double>("slot_Centroid_y").key;
    auto xsKey = Cat.getSchema().find<float>("slot_Centroid_xSigma").key;
    auto ysKey = Cat.getSchema().find<float>("slot_Centroid_ySigma").key;
    auto mxxKey = Cat.getSchema().find<double>("slot_Shape_xx").key;
    auto myyKey = Cat.getSchema().find<double>("slot_Shape_yy").key;
    auto mxyKey = Cat.getSchema().find<double>("slot_Shape_xy").key;
    auto fluxKey = Cat.getSchema().find<double>(fluxField + "_flux").key;
    auto efluxKey = Cat.getSchema().find<double>(fluxField  + "_fluxSigma").key;

    wholeCatalog.clear();
    for (auto i = Cat.begin(); i != Cat.end(); ++i)
    {
        MeasuredStar *ms = new MeasuredStar();
        ms->x = i->get(xKey);
        ms->y = i->get(yKey);
        ms->vx = sq(i->get(xsKey));
        ms->vy = sq(i->get(ysKey));
        /* the xy covariance is not provided in the input catalog: we
        cook it up from the x and y position variance and the shape
         measurements: */
        double mxx = i->get(mxxKey);
        double myy = i->get(myyKey);
        double mxy = i->get(mxyKey);
        ms->vxy = mxy*(ms->vx + ms->vy)/(mxx + myy);
        if (ms->vx < 0 || ms->vy < 0 || (ms->vxy*ms->vxy) > (ms->vx*ms->vy)) {
            std::cout << "Bad source detected in LoadCatalog : " << ms->vx << " " << ms->vy << " " <<
                      ms->vxy*ms->vxy << " " << ms->vx*ms->vy << std::endl;
            continue;
        }
        ms->flux = i->get(fluxKey);
        ms->eflux = i->get(efluxKey);
        ms->mag = -2.5*log10(ms->flux) + zp;
        ms->setCcdImage(this);
        wholeCatalog.push_back(ms);
    }
    wholeCatalog.setCcdImage(this);
}

CcdImage::CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
                   const Point &CommonTangentPoint,
                   const PTR(lsst::afw::image::TanWcs) wcs,
                   const PTR(lsst::daf::base::PropertySet) meta,
//                   const PTR(lsst::afw::image::VisitInfo) visitInfo,
                   const lsst::afw::geom::Box2I &bbox,
                   const std::string &filter,
                   const PTR(lsst::afw::image::Calib) calib,
                   const int &visit,
                   const int &ccdId,
                   const std::string &fluxField ) :

    _filter(filter), _visit(visit), _ccdId(ccdId), index(-1), expindex(-1),
    commonTangentPoint(CommonTangentPoint)

{
    // zero point
    zp = 2.5*log10(calib->getFluxMag0().first);

    LoadCatalog(Ri, fluxField);

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    imageFrame = Frame(lowerLeft, upperRight);

    readWcs = new jointcal::TanSipPix2RaDec(jointcal::ConvertTanWcs(wcs));

    // use some other variable in case we later have to actually convert the
    // read wcs:
    const BaseTanWcs* tanWcs = readWcs.get();

    inverseReadWcs = readWcs->InverseTransfo(0.01, imageFrame);

    std::stringstream out;
    out << visit << "_" << ccdId;
    riName = out.str();

    /* we don't assume here that we know the internals of TanPix2RaDec:
       to construct pix->TP, we do pix->sky->TP, although pix->sky
       actually goes through TP */

    GtransfoLin identity;
    TanRaDec2Pix raDec2TP(identity, tanWcs->TangentPoint());
    pix2TP = GtransfoCompose(&raDec2TP, tanWcs);
    TanPix2RaDec CTP2RaDec(identity, CommonTangentPoint);
    CTP2TP = GtransfoCompose(&raDec2TP, &CTP2RaDec);

    // jump from one TP to an other:
    TanRaDec2Pix raDec2CTP(identity, CommonTangentPoint);
    //  TanPix2RaDec TP2RaDec(identity, tanWcs->TangentPoint());
    //  TP2CTP = GtransfoCompose(&raDec2CTP, &TP2RaDec);
    TanPix2RaDec TP2RaDec(identity, tanWcs->TangentPoint());
    TP2CTP = GtransfoCompose(&raDec2CTP, &TP2RaDec);
    sky2TP = new TanRaDec2Pix(identity, tanWcs->TangentPoint());

    // this one is needed for matches :
    pix2CommonTangentPlane = GtransfoCompose(&raDec2CTP, tanWcs);

//    double latitude = visitInfo->getObservatory().getLatitude();
//    double lst_obs = visitInfo->getEra();
//    double ra = visitInfo->getBoresightRaDec().getRa();
//    double dec = visitInfo->getBoresightRaDec().getDec();
//    double hourAngle = visitInfo->getBoresightHourAngle();
//    airMass = visitInfo->getBoresightAirmass();
//    mjd = visitInfo->getDate().get(lsst::daf::base::DateTime::MJD);

//    airMass = meta->get<double>("AIRMASS");
//    mjd = meta->get<double>("MJD-OBS");  // Julian date
    airMass = meta->get<double>("BORE-AIRMASS");
    mjd = meta->get<double>("MJDATE");  // Julian date
//    expTime = meta->get<double>("EXPTIME");
//    double latitude = lsst::afw::geom::degToRad(meta->get<double>("LATITUDE"));
    double latitude = lsst::afw::geom::degToRad(meta->get<double>("OBS-LAT"));
//    double lst_obs = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("LST-OBS")));
//    double lst_obs = lsst::afw::geom::degToRad(meta->get<double>("AVG-ERA"));
//    double ra = lsst::afw::geom::degToRad(meta->get<double>("RA_DEG"));
//    double dec = lsst::afw::geom::degToRad(meta->get<double>("DEC_DEG"));
    double ra = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("RA")));
    double dec = lsst::afw::geom::degToRad(DecStringToDeg(meta->get<std::string>("DEC")));

    hourAngle = lsst::afw::geom::degToRad(meta->get<double>("AVG-ERA")+meta->get<double>("OBS-LONG")-meta->get<double>("BORE-RA"));

//    hourAngle = (lst_obs-ra);
//    std::cout << hourAngle << std::endl;
//    if  (hourAngle>M_PI) hourAngle -= 2*M_PI;
//    if  (hourAngle<-M_PI) hourAngle += 2*M_PI;

    // lsstSim doesn't manage ERA (and thus Hour Angle) properly, so it's going to be NaN.
    // Because we need the refraction vector later, go with 0 HA to prevent crashes on that NaN.
    if (std::isnan(hourAngle) == true) {
        hourAngle = 0;
    }

    if (airMass == 1)
        sineta = coseta = tgz = 0;
    else
    {
        double cosz = 1./airMass;
        double sinz = sqrt(1 - cosz*cosz); //astronomers usually observe above the horizon
        tgz = sinz/cosz;
        sineta = cos(latitude)*sin(hourAngle)/sinz;
        coseta = sqrt(1 - sineta*sineta);
        if (dec > latitude) coseta = -coseta;
    }
//        std::cout << airMass << " " << ra << " " << dec << " " << sineta << " " << hourAngle << " " << latitude << std::endl;
}

}
} // end of namespaces
