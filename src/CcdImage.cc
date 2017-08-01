#include <assert.h>
#include <string>
#include <sstream>
#include <math.h>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/Point.h"

#include "lsst/log/Log.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Point.h"

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.CcdImage");
}

namespace lsst {
namespace jointcal {

void CcdImage::LoadCatalog(afw::table::SourceCatalog const &catalog, std::string const &fluxField) {
    auto xKey = catalog.getSchema().find<double>("slot_Centroid_x").key;
    auto yKey = catalog.getSchema().find<double>("slot_Centroid_y").key;
    auto xsKey = catalog.getSchema().find<float>("slot_Centroid_xSigma").key;
    auto ysKey = catalog.getSchema().find<float>("slot_Centroid_ySigma").key;
    auto mxxKey = catalog.getSchema().find<double>("slot_Shape_xx").key;
    auto myyKey = catalog.getSchema().find<double>("slot_Shape_yy").key;
    auto mxyKey = catalog.getSchema().find<double>("slot_Shape_xy").key;
    auto fluxKey = catalog.getSchema().find<double>(fluxField + "_flux").key;
    auto fluxErrKey = catalog.getSchema().find<double>(fluxField + "_fluxSigma").key;

    _wholeCatalog.clear();
    for (auto const &record : catalog) {
        auto ms = std::make_shared<MeasuredStar>();
        ms->x = record.get(xKey);
        ms->y = record.get(yKey);
        ms->vx = std::pow(record.get(xsKey), 2);
        ms->vy = std::pow(record.get(ysKey), 2);
        /* the xy covariance is not provided in the input catalog: we
        cook it up from the x and y position variance and the shape
         measurements: */
        double mxx = record.get(mxxKey);
        double myy = record.get(myyKey);
        double mxy = record.get(mxyKey);
        ms->vxy = mxy * (ms->vx + ms->vy) / (mxx + myy);
        if (ms->vx < 0 || ms->vy < 0 || (ms->vxy * ms->vxy) > (ms->vx * ms->vy)) {
            LOGLS_WARN(_log, "Bad source detected in LoadCatalog : " << ms->vx << " " << ms->vy << " "
                                                                     << ms->vxy * ms->vxy << " "
                                                                     << ms->vx * ms->vy);
            continue;
        }
        ms->setId(record.getId());
        ms->setInstFlux(record.get(fluxKey));
        ms->setInstFluxErr(record.get(fluxErrKey));
        // TODO: the below lines will be less clumsy once DM-4044 is cleaned up and we can say:
        // TODO: instFluxToMaggies(ms->getInstFlux(), ms) (because ms will be derived from afw::geom::Point).
        afw::geom::Point<double, 2> point(ms->x, ms->y);
        auto flux = _photoCalib->instFluxToMaggies(ms->getInstFlux(), ms->getInstFluxErr(), point);
        ms->setFlux(flux.value);
        ms->setFluxErr(flux.err);
        ms->mag = _photoCalib->instFluxToMagnitude(ms->getInstFlux(), point);
        ms->setCcdImage(this);
        _wholeCatalog.push_back(std::move(ms));
    }
    _wholeCatalog.setCcdImage(this);
}

CcdImage::CcdImage(afw::table::SourceCatalog &catalog, std::shared_ptr<lsst::afw::image::TanWcs> wcs,
                   std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, afw::geom::Box2I const &bbox,
                   std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib, int visit,
                   int ccdId, std::string const &fluxField)
        : _ccdId(ccdId), _visit(visit), _photoCalib(photoCalib), _filter(filter) {
    LoadCatalog(catalog, fluxField);

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    _imageFrame = Frame(lowerLeft, upperRight);

    _readWcs.reset(new jointcal::TanSipPix2RaDec(jointcal::convertTanWcs(wcs)));
    _inverseReadWcs = _readWcs->inverseTransfo(0.01, _imageFrame);

    std::stringstream out;
    out << visit << "_" << ccdId;
    _name = out.str();

    _boresightRaDec = visitInfo->getBoresightRaDec();
    _airMass = visitInfo->getBoresightAirmass();
    _mjd = visitInfo->getDate().get(lsst::daf::base::DateTime::MJD);
    double latitude = visitInfo->getObservatory().getLatitude();
    _lstObs = visitInfo->getEra();
    _hourAngle = visitInfo->getBoresightHourAngle();

    // lsstSim doesn't manage ERA (and thus Hour Angle) properly, so it's going to be NaN.
    // Because we need the refraction vector later, go with 0 HA to prevent crashes on that NaN.
    if (std::isnan(_hourAngle) == true) {
        _hourAngle = 0;
    }

    if (_airMass == 1)
        _sineta = _coseta = _tgz = 0;
    else {
        double cosz = 1. / _airMass;
        double sinz = sqrt(1 - cosz * cosz);  // astronomers usually observe above the horizon
        _tgz = sinz / cosz;
        _sineta = cos(latitude) * sin(_hourAngle) / sinz;
        _coseta = sqrt(1 - _sineta * _sineta);
        if (_boresightRaDec.getDec() > latitude) _coseta = -_coseta;
    }
}

void CcdImage::setCommonTangentPoint(Point const &commonTangentPoint) {
    _commonTangentPoint = commonTangentPoint;

    // use some other variable in case we later have to actually convert the
    // as-read wcs:
    const BaseTanWcs *tanWcs = _readWcs.get();

    /* we don't assume here that we know the internals of TanPix2RaDec:
       to construct pix->TP, we do pix->sky->TP, although pix->sky
       actually goes through TP */
    GtransfoLin identity;
    TanRaDec2Pix raDec2TP(identity, tanWcs->getTangentPoint());
    _pix2TP = gtransfoCompose(&raDec2TP, tanWcs);
    TanPix2RaDec CTP2RaDec(identity, commonTangentPoint);
    _CTP2TP = gtransfoCompose(&raDec2TP, &CTP2RaDec);

    // jump from one TP to an other:
    TanRaDec2Pix raDec2CTP(identity, commonTangentPoint);
    TanPix2RaDec TP2RaDec(identity, tanWcs->getTangentPoint());
    _TP2CTP = gtransfoCompose(&raDec2CTP, &TP2RaDec);
    _sky2TP.reset(new TanRaDec2Pix(identity, tanWcs->getTangentPoint()));

    // this one is needed for matches :
    _pix2CommonTangentPlane = gtransfoCompose(&raDec2CTP, tanWcs);
}
}  // namespace jointcal
}  // namespace lsst
