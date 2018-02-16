#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE test_trans

//The boost unit test header
#include "boost/test/unit_test.hpp"

#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/image/Image.h"
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
  auto propSet = lsst::afw::fits::readMetadata(fileName);
  auto skyWcs = lsst::afw::geom::makeSkyWcs(*propSet);
  jointcal::GtransfoSkyWcs gtransfoWcs(skyWcs);

  jointcal::Point where(100.,200.);
  jointcal::Point outPol = gtransfoWcs.apply(where);
  std::cout << std::setprecision(12) << "Poloka : " << outPol.x << ' ' << outPol.y << std::endl;

  lsst::afw::geom::Point2D whereSame(100.,200.);
  auto skyPos = skyWcs->pixelToSky(whereSame);
  lsst::afw::geom::Point2D outDeg = skyPos.getPosition(lsst::afw::geom::degrees);
  std::cout << "Stack : " << outDeg[0] << ' ' << outDeg[1] << std::endl;

  BOOST_CHECK_CLOSE(outPol.x, outDeg[0], .000001);
  BOOST_CHECK_CLOSE(outPol.y, outDeg[1], .000001);
}



/* test the GtransfoPoly::fit routine */

BOOST_AUTO_TEST_CASE(test_polyfit)
{
  std::string fileName = "tests/header_only.fits";

  lsst::afw::fits::Fits file(fileName, "r",0);
  auto propSet = lsst::afw::fits::readMetadata(fileName);
  auto skyWcs = lsst::afw::geom::makeSkyWcs(*propSet);
  jointcal::GtransfoSkyWcs gtransfoWcs(skyWcs);

  jointcal::StarMatchList sml;
  jointcal::BaseStarList bsl1, bsl2;

  for (double x=10; x<2000; x+= 120)
    for (double y=20; y<4000; y+=160)
      {
	auto s1 = std::make_shared<jointcal::BaseStar>(x,y,1,0.01);
	s1->vx = 0.1;
	s1->vy = 0.2;
	s1->vxy = 0.05;
	auto s2 = std::make_shared<jointcal::BaseStar>();
	gtransfoWcs.transformPosAndErrors(*s1, *s2);
	bsl1.push_back(s1);
	bsl2.push_back(s2);
	sml.push_back(jointcal::StarMatch(*s1,*s2,s1,s2));
      }
  jointcal::GtransfoPoly pol(3);
  double chi2 = pol.fit(sml);
  std::cout << " chi2/ndf " << chi2 << '/' << sml.size()-pol.getNpar() << std::endl;
  // since there is no noise, the chi2 should be very very small:
  BOOST_CHECK( fabs(chi2)<1e-8);
}


BOOST_AUTO_TEST_SUITE_END()
