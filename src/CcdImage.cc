#include <assert.h>
#include <string>
#include <sstream>
#include <math.h>

#include "lsst/log/Log.h"
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

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.CcdImage");
}

namespace lsst {
namespace jointcal {

static double sq(double x) { return x*x;}

void CcdImage::LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &catalog,
                           const std::string &fluxField)
{
    auto xKey = catalog.getSchema().find<double>("slot_Centroid_x").key;
    auto yKey = catalog.getSchema().find<double>("slot_Centroid_y").key;
    auto xsKey = catalog.getSchema().find<float>("slot_Centroid_xSigma").key;
    auto ysKey = catalog.getSchema().find<float>("slot_Centroid_ySigma").key;
    auto mxxKey = catalog.getSchema().find<double>("slot_Shape_xx").key;
    auto myyKey = catalog.getSchema().find<double>("slot_Shape_yy").key;
    auto mxyKey = catalog.getSchema().find<double>("slot_Shape_xy").key;
    auto fluxKey = catalog.getSchema().find<double>(fluxField + "_flux").key;
    auto efluxKey = catalog.getSchema().find<double>(fluxField  + "_fluxSigma").key;

    _wholeCatalog.clear();
    for (auto const &i: catalog)
    {
        MeasuredStar *ms = new MeasuredStar();
        ms->x = i.get(xKey);
        ms->y = i.get(yKey);
        ms->vx = sq(i.get(xsKey));
        ms->vy = sq(i.get(ysKey));
        /* the xy covariance is not provided in the input catalog: we
        cook it up from the x and y position variance and the shape
         measurements: */
        double mxx = i.get(mxxKey);
        double myy = i.get(myyKey);
        double mxy = i.get(mxyKey);
        ms->vxy = mxy*(ms->vx + ms->vy)/(mxx + myy);
        if (ms->vx < 0 || ms->vy < 0 || (ms->vxy*ms->vxy) > (ms->vx*ms->vy)) {
            LOGLS_WARN(_log, "Bad source detected in LoadCatalog : " << ms->vx << " " << ms->vy << " "
                       << ms->vxy*ms->vxy << " " << ms->vx*ms->vy);
            continue;
        }
        ms->flux = i.get(fluxKey);
        ms->eflux = i.get(efluxKey);
        ms->mag = _calib->getMagnitude(ms->flux);
        ms->setCcdImage(this);
        _wholeCatalog.push_back(ms);
    }
    _wholeCatalog.setCcdImage(this);
}

CcdImage::CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
                   const PTR(lsst::afw::image::TanWcs) wcs,
                   const PTR(lsst::afw::image::VisitInfo) visitInfo,
                   const lsst::afw::geom::Box2I &bbox,
                   const std::string &filter,
                   const PTR(lsst::afw::image::Calib) calib,
                   const int &visit,
                   const int &ccdId,
                   const std::string &fluxField ) :

    _filter(filter), _visit(visit), _ccdId(ccdId), _calib(calib)

{
    LoadCatalog(Ri, fluxField);

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    imageFrame = Frame(lowerLeft, upperRight);

    readWcs = new jointcal::TanSipPix2RaDec(jointcal::convertTanWcs(wcs));
    inverseReadWcs = readWcs->InverseTransfo(0.01, imageFrame);

    std::stringstream out;
    out << visit << "_" << ccdId;
    name = out.str();

    boresightRaDec = visitInfo->getBoresightRaDec();
    airMass = visitInfo->getBoresightAirmass();
    mjd = visitInfo->getDate().get(lsst::daf::base::DateTime::MJD);
    double latitude = visitInfo->getObservatory().getLatitude();
    double lst_obs = visitInfo->getEra();
    double hourAngle = visitInfo->getBoresightHourAngle();

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
        if (boresightRaDec.getDec() > latitude) coseta = -coseta;
    }
}

void CcdImage::setCommonTangentPoint(const Point &commonTangentPoint)
{
    _commonTangentPoint = commonTangentPoint;

    // use some other variable in case we later have to actually convert the
    // as-read wcs:
    const BaseTanWcs* tanWcs = readWcs.get();

    /* we don't assume here that we know the internals of TanPix2RaDec:
       to construct pix->TP, we do pix->sky->TP, although pix->sky
       actually goes through TP */
    GtransfoLin identity;
    TanRaDec2Pix raDec2TP(identity, tanWcs->TangentPoint());
    pix2TP = GtransfoCompose(&raDec2TP, tanWcs);
    TanPix2RaDec CTP2RaDec(identity, commonTangentPoint);
    CTP2TP = GtransfoCompose(&raDec2TP, &CTP2RaDec);

    // jump from one TP to an other:
    TanRaDec2Pix raDec2CTP(identity, commonTangentPoint);
    TanPix2RaDec TP2RaDec(identity, tanWcs->TangentPoint());
    TP2CTP = GtransfoCompose(&raDec2CTP, &TP2RaDec);
    sky2TP = new TanRaDec2Pix(identity, tanWcs->TangentPoint());

    // this one is needed for matches :
    pix2CommonTangentPlane = GtransfoCompose(&raDec2CTP, tanWcs);
}

}
} // end of namespaces
