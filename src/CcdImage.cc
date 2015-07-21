#include <assert.h>
#include <string>

#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/SipToGtransfo.h"
#include "lsst/meas/simastrom/AstroUtils.h"
#include "lsst/afw/image/Image.h"
#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/meas/simastrom/Point.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Angle.h"

namespace simAstrom = lsst::meas::simastrom;
namespace afwImg = lsst::afw::image;

namespace lsst {
namespace meas {
namespace simastrom {

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

void CcdImage::LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat)
{
  auto xKey = Cat.getSchema().find<double>("base_SdssCentroid_x").key;
  auto yKey = Cat.getSchema().find<double>("base_SdssCentroid_y").key;
  auto xsKey = Cat.getSchema().find<float>("base_SdssCentroid_xSigma").key;
  auto ysKey = Cat.getSchema().find<float>("base_SdssCentroid_ySigma").key;
  auto mxxKey = Cat.getSchema().find<double>("base_SdssShape_xx").key;
  auto myyKey = Cat.getSchema().find<double>("base_SdssShape_yy").key;
  auto mxyKey = Cat.getSchema().find<double>("base_SdssShape_xy").key;
  auto fluxKey = Cat.getSchema().find<double>("base_CircularApertureFlux_5_flux").key;
  auto efluxKey = Cat.getSchema().find<double>("base_CircularApertureFlux_5_fluxSigma").key;
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
      ms->flux = i->get(fluxKey);
      ms->eflux = i->get(efluxKey);
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
            const PTR(lsst::afw::image::Calib) calib  ) :
            
    index(-1), expindex(-1),
    commonTangentPoint(CommonTangentPoint)

{
  LoadCatalog(Ri);
      // Just checking that we get something sensible
//    std::cout << Ri[10].getRa() << std::endl;
//    std::cout << wcs->getPixelOrigin() << std::endl;
//    std::cout << meta->get<double>("LATITUDE") << std::endl;
//    std::cout << bbox << std::endl;
//    std::cout << filter << std::endl;

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    imageFrame = Frame(lowerLeft, upperRight);
    
    readWcs = new simAstrom::TanSipPix2RaDec(simAstrom::ConvertTanWcs(wcs));
    
    // Don't know if the following is necessary. Why shouldn't we use readWcs everywhere instead of tanWcs ?
    TanSipPix2RaDec *tanWcs = dynamic_cast<TanSipPix2RaDec*>(readWcs);
  
    inverseReadWcs = readWcs->InverseTransfo(0.01, imageFrame);
 
    band = filter;
    bandIndex = getBandIndex(band);

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
    
    // In the following we read informations directly from the fits header. This will have to be modified in order to be
    // instrument independent
    
    airMass = meta->get<double>("AIRMASS");
    jd = meta->get<double>("MJD-OBS");  // Julian date
    expTime = meta->get<double>("EXPTIME");
    double latitude = lsst::afw::geom::degToRad(meta->get<double>("LATITUDE"));
    double lst_obs = lsst::afw::geom::degToRad(RaStringToDeg(meta->get<std::string>("LST-OBS")));
    double ra = lsst::afw::geom::degToRad(meta->get<double>("RA_DEG"));
    double dec = lsst::afw::geom::degToRad(meta->get<double>("DEC_DEG"));
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
    
}}} // end of namespaces
