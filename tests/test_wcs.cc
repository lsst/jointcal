#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_wcs
 
//The boost unit test header
#include "boost/test/unit_test.hpp"

#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/meas/simastrom/SipToGtransfo.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/fits.h"
#include "lsst/daf/base.h" 

namespace simAstrom = lsst::meas::simastrom; 
namespace afwImg = lsst::afw::image;

/* Test TanSipPix2RaDec against pixelToSky */

BOOST_AUTO_TEST_CASE(test_wcs)
{
    
  std::string fileName = "tests/header_only.fits";

  lsst::afw::fits::Fits file(fileName, "r",0);
  PTR(lsst::daf::base::PropertySet) propSet = afwImg::readMetadata(fileName);
  PTR(afwImg::Wcs) wcs = afwImg::makeWcs(propSet);
  
  const PTR(afwImg::TanWcs) tanWcs = boost::dynamic_pointer_cast<afwImg::TanWcs>(wcs);
  
  simAstrom::TanSipPix2RaDec gtransfoWcs = simAstrom::ConvertTanWcs(tanWcs);
  simAstrom::Point where(100.,200.);
  simAstrom::Point outPol = gtransfoWcs.apply(where);
  std::cout << std::setprecision(12) << "Poloka : " << outPol.x << ' ' << outPol.y << std::endl;
  
  lsst::afw::geom::Point2D whereSame(100.,200.);
  PTR(lsst::afw::coord::Coord) coord = wcs->pixelToSky(whereSame);
  lsst::afw::geom::Point2D outDeg = coord->getPosition(lsst::afw::geom::degrees);
  std::cout << "Stack : " << outDeg[0] << ' ' << outDeg[1] << std::endl;
  
  BOOST_CHECK_CLOSE(outPol.x, outDeg[0], .000001);
  BOOST_CHECK_CLOSE(outPol.y, outDeg[1], .000001);

}
