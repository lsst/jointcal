#include <assert.h>
#include <string>

#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/SipToGtransfo.h"
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
    
double RaStringToDeg(const std::string RaString)
{
int hours, minutes; double seconds;
if (sscanf(RaString.c_str(),"%d:%d:%lf", &hours, &minutes, &seconds) == 3)
  {
    return 15.*( double(hours) + double(minutes)/60. + seconds/3600.);
  }
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
    commonTangentPoint(CommonTangentPoint),
    wholeCatalog(Ri)

{
      // Just checking that we get something sensible
    std::cout << Ri[10].getRa() << std::endl;
    std::cout << wcs->getPixelOrigin() << std::endl;
    std::cout << meta->get<double>("LATITUDE") << std::endl;
    std::cout << bbox << std::endl;
    std::cout << filter << std::endl;

    Point lowerLeft(bbox.getMinX(), bbox.getMinY());
    Point upperRight(bbox.getMaxX(), bbox.getMaxY());
    imageFrame = Frame(lowerLeft, upperRight);
    
    simAstrom::TanSipPix2RaDec tt = simAstrom::ConvertTanWcs(wcs);
    readWcs = &tt;
  
    inverseReadWcs = readWcs->InverseTransfo(0.01, imageFrame);
    
    band = filter;
    bandIndex = getBandIndex(band);
    
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
}
    
}}} // end of namespaces
