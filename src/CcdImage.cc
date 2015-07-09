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
    
double RaStringToDeg(const string RaString)
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
            const std::string &filter ) :
            
    index(-1), expindex(-1),
    commonTangentPoint(CommonTangentPoint)
  
  
{
      // Just checking that we get something sensible
    std::cout << Ri[10].getRa() << std::endl;
    std::cout << wcs->getPixelOrigin() << std::endl;
    std::cout << meta->get<double>("LATITUDE") << std::endl;
    std::cout << bbox << std::endl;
    std::cout << filter << std::endl;
    
  // needed transfos
  // FitsHeader head(Ri.FitsName());

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
    double ra = lsst::afw::geom::degToRad(meta->get<double>("RA_DEG"));
    double dec = lsst::afw::geom::degToRad(meta->get<double>("DEC_DEG"));
//    double lst_obs = RaStringToDeg(meta->get
}
  
#ifdef TO_BE_FIXED   

  chip = head.KeyVal("TOADCHIP");
  band = string(head.KeyVal("TOADBAND"));
  instrument = string(head.KeyVal("TOADINST"));
  bandRank = 0; // will be set by Associations if pertinent. 
  bandIndex = getBandIndex(band);
  if (head.HasKey("EXPNUM"))
    shoot = head.KeyVal("EXPNUM"); // specific to Megacam
  else if (head.HasKey("IDPROCES"))
    { // UHSnifs
      string idproc = head.KeyVal("IDPROCESS");
      sscanf(idproc.c_str(),"%11u",&shoot);
    }
      
  jd = Ri.ModifiedJulianDate();

  gfseeing = Ri.GFSeeing();
  seeing = Ri.Seeing();
  sigmaback = Ri.SigmaBack();

  flatName = string(head.KeyVal("IMRED_FF")); // e.g. 04Am02.flat.i.36.02.fits[ccd06]
  flatName = flatName.substr(0,flatName.find('[')); // remove [ccdxx]
  
  // update the real flats
  char buffer[PATH_MAX];
  string flatpath = Preferences().flatpath;
  //  flatpath = flatpath + "/flat/" + band + "/" + flatName.substr(0,flatName.find(".fits"));
  flatpath = flatpath + "/flat/" + flatName.substr(0,flatName.find(".fits"));
  //  sprintf(buffer, 
  //	  "%s/flat/%s/%s/ccd_%02d.fits", 
  //	  flatpath.c_str(), band.c_str(), flatName.substr(0,flatName.find(".fits"));
  
  // fetch the elixir scatter map 
  //  cfhtscatter = flatpath + "/scatterflat";
  cfhtscatter = flatpath + "/scatter";
  char* ret = realpath(cfhtscatter.c_str(), buffer);
  if(ret)
    cfhtscatter = buffer;
  else
    cfhtscatter = "";
  
  // get the SNLS grid corrections 
  // -- same as the elixir corrections, 
  //    but determined here. 
  //  snlsgrid = flatpath + "/snlsgrid";
  snlsgrid = flatpath + "/snlsgrid-V0";
  ret = realpath(snlsgrid.c_str(), buffer);
  if(ret)
    snlsgrid = buffer;
  else
    snlsgrid = "";
  
  // get the flatcvmap 
  flatcvmap = flatpath + "/cvmap";
  ret = realpath(flatcvmap.c_str(), buffer);
  if(ret)
    flatcvmap = buffer;
  else
    flatcvmap = "";
  
  sprintf(buffer, "%s/mosaic_4_9_%02d.fits", cfhtscatter.c_str(), chip);
  cfhtscatter = buffer;
  sprintf(buffer, "%s/mosaic_4_9_%02d.fits", snlsgrid.c_str(), chip);
  snlsgrid = buffer;
  sprintf(buffer, "%s/mosaic_4_9_%02d.fits", flatcvmap.c_str(), chip);
  flatcvmap = buffer;
  
  if(!FileExists(cfhtscatter)) cfhtscatter="";
  if(!FileExists(snlsgrid)) snlsgrid="";
  if(!FileExists(flatcvmap)) flatcvmap="";
  
//   cout << " flatpath=" << flatpath << endl
//        << " cfhtscatter=" << cfhtscatter << endl
//        << " snlsgrid=" << snlsgrid << endl
//        << " flatcvmap=" << flatcvmap << endl;
  
  airMass = (double)head.KeyVal("TOADAIRM");
  //  dateObs = head.KeyVal("TOADDATE");
  dateObs = (string)head.KeyVal("DATE-OBS"); // to have the same DATE format as in the DATABASE
  expTime = head.KeyVal("TOADEXPO");

  /* Megacam fits header describing how to use Elixir photm info:

COMMENT   Formula for Photometry, based on keywords given in this header:
COMMENT   m = -2.5*log(DN) + 2.5*log(EXPTIME)
COMMENT   M = m + PHOT_C + PHOT_K*(AIRMASS - 1) + PHOT_X*(PHOT_C1 - PHOT_C2)

we don't use the color terms since we are working in a single band.
  */

  photk = head.KeyVal("PHOT_K");
  fluxCoeff = pow(10., -0.4*(photk*(airMass-1)))/expTime;
  photc = head.KeyVal("PHOT_C");
  elixirZP = photc/*+photk*(airMass-1)+2.5*log10(expTime)*/;
  zp = head.KeyVal(ZpKey);
  DataCards psfzpFile(Ri.Dir()+"/psfzp.dat");
  if(psfzpFile.HasKey("PSFZP")) psfzp = psfzpFile.DParam("PSFZP");
  else psfzp = zp;

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

  // pick up what is needed for atmospheric refraction
  double airmass = head.KeyVal("AIRMASS");
  double latitude = 19.825252*M_PI/180;
  if (head.HasKey("LATITUDE"))
    latitude= double(head.KeyVal("LATITUDE"))*M_PI/180;
  else
    {
      cout << " no latitude fits key for image " << Ri.Name() << " hope it was taken at CFHT" << endl;
    }
      
  double lst_obs= RaStringToDeg(head.KeyVal("LST-OBS"))*M_PI/180;
  double ra = RaStringToDeg(head.KeyVal("RA"))*M_PI/180;
  double dec = DecStringToDeg(head.KeyVal("DEC"))*M_PI/180;
  hourAngle = (lst_obs-ra);
  if  (hourAngle>M_PI) hourAngle -= 2*M_PI;
  if  (hourAngle<-M_PI) hourAngle += 2*M_PI;

  if (airmass==1)
    sineta = coseta = tgz = 0;

  else
    {
      double cosz = 1./airmass;
      double sinz = sqrt(1-sqr(cosz)); //astronomers usually observe above the horizon 
      tgz = sinz/cosz;
      sineta = cos(latitude)*sin(hourAngle)/sinz;
      coseta = sqrt(1-sqr(sineta));
      if (dec > latitude) coseta = -coseta;
    }


  toadsZeroPoint = Ri.ZeroPoint();

  // read the catalog
  LoadIt->load(Ri, *this, wholeCatalog);
  if (wholeCatalog.size() == 0)
    {
      cerr << " empty catalog for " << Ri.Name() << endl;
      return;
    }
  wholeCatalog.SetCcdImage(this);
}
#endif
    
    
}}} // end of namespaces
