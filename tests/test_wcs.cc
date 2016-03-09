#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE test_trans
 
//The boost unit test header
#include "boost/test/unit_test.hpp"

#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/fits.h"
#include "lsst/daf/base.h" 

#include <stdlib.h> /* for getenv */

// NOTE: turn this flag on to raise exceptions on floating point errors.
// NOTE: this only works on GNU/Linux (fenableexcept is not C++ standard).
// #define DUMP_CORE_ON_FPE
#ifdef DUMP_CORE_ON_FPE
#define _GNU_SOURCE 1
#define __USE_GNU
#include <fenv.h>
static void __attribute__ ((constructor)) trapfpe ()
{
   // Enable some exceptions.  At startup all exceptions are masked.
  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif


namespace jointcal = lsst::jointcal; 
namespace afwImg = lsst::afw::image;

/* Test jointcal::TanSipPix2RaDec::apply against afwImg::Wcs::pixelToSky */

BOOST_AUTO_TEST_SUITE(test_transfos)

BOOST_AUTO_TEST_CASE(test_wcs)
{
  
  std::string fileName = "tests/header_only.fits";

  lsst::afw::fits::Fits file(fileName, "r",0);
  PTR(lsst::daf::base::PropertySet) propSet = afwImg::readMetadata(fileName);
  PTR(afwImg::Wcs) wcs = afwImg::makeWcs(propSet);
  
  const PTR(afwImg::TanWcs) tanWcs = boost::dynamic_pointer_cast<afwImg::TanWcs>(wcs);
  
  jointcal::TanSipPix2RaDec gtransfoWcs = jointcal::ConvertTanWcs(tanWcs);
  jointcal::Point where(100.,200.);
  jointcal::Point outPol = gtransfoWcs.apply(where);
  std::cout << std::setprecision(12) << "Poloka : " << outPol.x << ' ' << outPol.y << std::endl;
  
  lsst::afw::geom::Point2D whereSame(100.,200.);
  PTR(lsst::afw::coord::Coord) coord = wcs->pixelToSky(whereSame);
  lsst::afw::geom::Point2D outDeg = coord->getPosition(lsst::afw::geom::degrees);
  std::cout << "Stack : " << outDeg[0] << ' ' << outDeg[1] << std::endl;
  
  BOOST_CHECK_CLOSE(outPol.x, outDeg[0], .000001);
  BOOST_CHECK_CLOSE(outPol.y, outDeg[1], .000001);

}



/* test the GtransfoPoly::fit routine */

BOOST_AUTO_TEST_CASE(test_polyfit)
{
  std::string fileName = "tests/header_only.fits";

  lsst::afw::fits::Fits file(fileName, "r",0);
  PTR(lsst::daf::base::PropertySet) propSet = afwImg::readMetadata(fileName);
  PTR(afwImg::Wcs) wcs = afwImg::makeWcs(propSet);
  
  const PTR(afwImg::TanWcs) tanWcs = boost::dynamic_pointer_cast<afwImg::TanWcs>(wcs);
  
  jointcal::TanSipPix2RaDec gtransfoWcs = jointcal::ConvertTanWcs(tanWcs);

  jointcal::StarMatchList sml;
  jointcal::BaseStarList bsl1, bsl2;

  for (double x=10; x<2000; x+= 120)
    for (double y=20; y<4000; y+=160)
      {
	jointcal::BaseStar *s1 = new jointcal::BaseStar(x,y,1);
	s1->vx = 0.1;
	s1->vy = 0.2;
	s1->vxy = 0.05;
	jointcal::BaseStar *s2 = new jointcal::BaseStar();
	gtransfoWcs.TransformPosAndErrors(*s1, *s2);
	bsl1.push_back(s1);
	bsl2.push_back(s2);		       
	sml.push_back(jointcal::StarMatch(*s1,*s2,s1,s2));
      }
  jointcal::GtransfoPoly pol(3);
  double chi2 = pol.fit(sml);
  std::cout << " chi2/ndf " << chi2 << '/' << sml.size()-pol.Npar() << std::endl;
  // since there is no noise, the chi2 should be very very small:
  BOOST_CHECK( fabs(chi2)<1e-8);
}


/* this routine checks that converting a WCS from afs to jointcal and
   back to afw does not change anything. It checks both that pix2sky and 
   sky2pix transformations are preserved. */
BOOST_AUTO_TEST_CASE(test_wcs_convertions)
{
  std::string fileName = "tests/header_only.fits";

  lsst::afw::fits::Fits file(fileName, "r",0);
  PTR(lsst::daf::base::PropertySet) propSet = afwImg::readMetadata(fileName);
  PTR(afwImg::Wcs) wcs = afwImg::makeWcs(propSet);
  
  const PTR(afwImg::TanWcs) tanWcs = boost::dynamic_pointer_cast<afwImg::TanWcs>(wcs);
  
  jointcal::TanSipPix2RaDec gtransfoWcs = jointcal::ConvertTanWcs(tanWcs);
  int naxis1 = propSet->get<int>("NAXIS1");
  int naxis2 = propSet->get<int>("NAXIS2");
  jointcal::Frame imageFrame(0,0,naxis1,naxis2);



  // test the back conversion, in two cases
  for (int noLowOrderSipTerm  = 0; noLowOrderSipTerm <=1; noLowOrderSipTerm++) {
    PTR(afwImg::TanWcs) tanWcs2 = GtransfoToTanWcs(gtransfoWcs, 
						   imageFrame,
						   noLowOrderSipTerm);
    lsst::afw::geom::Point2D where(1000.,200.);

    PTR(lsst::afw::coord::Coord) coord = wcs->pixelToSky(where);
    /* one afw check : if we use pix->sky->pix with an afw WCS, how
       close do we get back ?  */

    lsst::afw::geom::Point2D whereBack = wcs->skyToPixel((*coord)[0], (*coord)[1]);
    BOOST_CHECK_CLOSE(where[0], whereBack[0], 1e-3);
    BOOST_CHECK_CLOSE(where[1], whereBack[1], 1e-3);


    lsst::afw::geom::Point2D outDeg = coord->getPosition(lsst::afw::geom::degrees);

    PTR(lsst::afw::coord::Coord) coord2 = tanWcs2->pixelToSky(where);
    lsst::afw::geom::Point2D outDeg2 = coord2->getPosition(lsst::afw::geom::degrees);
    // Check the proximity on sky
    BOOST_CHECK_CLOSE(outDeg[0], outDeg2[0], 1e-7);
    BOOST_CHECK_CLOSE(outDeg[1], outDeg2[1], 1e-7);

    // sky to pixels
    lsst::afw::geom::Point2D whereBack2 = tanWcs2->skyToPixel((*coord2)[0], (*coord2)[1]);


    BOOST_CHECK_CLOSE(whereBack[0], whereBack2[0], 1e-4);
    BOOST_CHECK_CLOSE(whereBack[1], whereBack2[1], 1e-4);
  }


}


BOOST_AUTO_TEST_SUITE_END()
