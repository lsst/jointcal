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

  /* Interesting fields of the stack catalogs :
 'base_SdssCentroid_x'
 'base_SdssCentroid_y'
 'base_SdssCentroid_xSigma'
 'base_SdssCentroid_ySigma'

  We miss the xy uncertainty term.
  We can cook it up from the sdss shape:
  'base_SdssShape_xx'
  'base_SdssShape_yy'
  'base_SdssShape_xy'

  for fluxes, we might use :
  'base_CircularApertureFlux_2_flux'
  'base_CircularApertureFlux_2_fluxSigma'

   where the '2' should be read from the environment.
  */

static double sq(const double &x) { return x*x;}

void CcdImage::LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat,const std::string &fluxField)
{
  auto xKey = Cat.getSchema().find<double>("base_SdssCentroid_x").key;
  auto yKey = Cat.getSchema().find<double>("base_SdssCentroid_y").key;
  auto xsKey = Cat.getSchema().find<float>("base_SdssCentroid_xSigma").key;
  auto ysKey = Cat.getSchema().find<float>("base_SdssCentroid_ySigma").key;
  auto mxxKey = Cat.getSchema().find<double>("base_SdssShape_xx").key;
  auto myyKey = Cat.getSchema().find<double>("base_SdssShape_yy").key;
  auto mxyKey = Cat.getSchema().find<double>("base_SdssShape_xy").key;
  auto fluxKey = Cat.getSchema().find<double>(fluxField + "_flux").key;
  auto efluxKey = Cat.getSchema().find<double>(fluxField  + "_fluxSigma").key;
  
  wholeCatalog.clear();
  for (auto i = Cat.begin(); i !=Cat.end(); ++i)
    {
      MeasuredStar *ms = new MeasuredStar();
      ms->x = i->get(xKey);
      ms->y = i->get(yKey);
      ms->vx = sq(i->get(xsKey));
      ms->vy = sq(i->get(ysKey));
      /* the xy covariance is not provided in the input catalog: we
	 cook it up from the x and y position variance and the shape
	 measurements: */
      double mxx= i->get(mxxKey);
      double myy= i->get(myyKey);
      double mxy= i->get(mxyKey);
      ms->vxy = mxy*(ms->vx+ms->vy)/(mxx+myy);
      if (ms->vx < 0 || ms->vy< 0 || (ms->vxy*ms->vxy)>(ms->vx*ms->vy)) {
          std::cout << "Bad source detected in LoadCatalog : " << ms->vx << " " << ms->vy << " " <<
          ms->vxy*ms->vxy << " " << ms->vx*ms->vy << std::endl;
          continue;
        }
      ms->flux = i->get(fluxKey);
      ms->eflux = i->get(efluxKey);
      ms->mag = -2.5*log10(ms->flux) + zp;
      ms->SetCcdImage(this);
      wholeCatalog.push_back(ms);
    }
  wholeCatalog.SetCcdImage(this);
}


    
static int getBandIndex(std::string const& band)
{
  if(band == "u") return 0;
  if(band == "g") return 1;
  if(band == "r") return 2;
  if(band == "i") return 3;
  if(band == "z") return 4;

  if(band == "U") return 5;
  if(band == "B") return 6;
  if(band == "V") return 7;
  if(band == "R") return 8;
  if(band == "I") return 9;
  return -1;
}


CcdImage::CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
            const Point &CommonTangentPoint,
            const PTR(lsst::afw::image::TanWcs) wcs,
            const PTR(lsst::daf::base::PropertySet) meta,
            const lsst::afw::geom::Box2I &bbox,
            const std::string &filter,
            const PTR(lsst::afw::image::Calib) calib,
            const int &visit,
            const int &ccd,
            const std::string &camera,
            const std::string &fluxField ) :
            
    index(-1), expindex(-1),
    commonTangentPoint(CommonTangentPoint)

{
    // zero point
    zp = 2.5*log10(calib->getFluxMag0().first);
    
    LoadCatalog(Ri, fluxField);
      // Just checking that we get something sensible
//    std::cout << Ri[10].getRa() << std::endl;
//    std::cout << wcs->getPixelOrigin() << std::endl;
//    std::cout << meta->get<double>("LATITUDE") << std::endl;
//    std::cout << bbox << std::endl;
//    std::cout << filter << std::endl;

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    imageFrame = Frame(lowerLeft, upperRight);
    
    readWcs = new jointcal::TanSipPix2RaDec(jointcal::ConvertTanWcs(wcs));

    // use some other variable in case we later have to actually convert the
    // read wcs:
    const BaseTanWcs* tanWcs = readWcs.get();

    inverseReadWcs = readWcs->InverseTransfo(0.01, imageFrame);
 
    band = filter;
    bandIndex = getBandIndex(band);
    chip = ccd;
    shoot = visit;
    
    std::stringstream out;
    out << visit << "_" << ccd;
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
    
    // In the following we read informations directly from the fits header which is instrument dependent
    // We rely on the camera name which is not optimal as a camera can be mounted on different telescopes.
    // We would rather need the telescope name.
    // This instrument specific part will eventually have to be handled by the camera mapper
    
    double latitude;
    double lst_obs;
    double ra;
    double dec;
    if (camera == "megacam") {
        airMass = meta->get<double>("AIRMASS");
        jd = meta->get<double>("MJD-OBS");  // Julian date
        expTime = meta->get<double>("EXPTIME");
        latitude = lsst::afw::geom::degToRad(meta->get<double>("LATITUDE"));
        lst_obs = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("LST-OBS")));
        ra = lsst::afw::geom::degToRad(meta->get<double>("RA_DEG"));
        dec = lsst::afw::geom::degToRad(meta->get<double>("DEC_DEG"));
    }
    else if (camera == "hsc") {
        airMass = meta->get<double>("AIRMASS");
        jd = meta->get<double>("MJD");  // Julian date
        expTime = meta->get<double>("EXPTIME");
        latitude = lsst::afw::geom::degToRad(19.825252);   // Does not seem to be in the header
        lst_obs = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("LST-STR")));
        ra = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("RA2000")));
        dec = lsst::afw::geom::degToRad(DecStringToDeg(meta->get<std::string>("DEC2000")));
    }
    else {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"CcdImage, unsupported camera "+camera);
    }
    hourAngle = (lst_obs-ra);
    if  (hourAngle>M_PI) hourAngle -= 2*M_PI;
    if  (hourAngle<-M_PI) hourAngle += 2*M_PI;
    
    if (airMass==1)
       sineta = coseta = tgz = 0;
    else
    {
      double cosz = 1./airMass;
      double sinz = sqrt(1-cosz*cosz); //astronomers usually observe above the horizon
      tgz = sinz/cosz;
      sineta = cos(latitude)*sin(hourAngle)/sinz;
      coseta = sqrt(1-sineta*sineta);
      if (dec > latitude) coseta = -coseta;
    }
    bandRank = 0; // will be set by Associations if pertinent.
}
    
}} // end of namespaces
