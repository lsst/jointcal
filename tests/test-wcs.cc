#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/fits.h"
#include "lsst/daf/base.h" 
namespace simAstrom = lsst::meas::simastrom; 
namespace afwImg = lsst::afw::image;

simAstrom::TanSipPix2RaDec ConvertTanWcs(const lsst::afw::image::TanWcs* wcs);

int main(int nargs, char **args)
{
  if (nargs<2)
    {
      std::cout << "need a file name " << std::endl;
      exit(1);
    }

  std::string fileName = args[1];

  lsst::afw::fits::Fits file(fileName, "r",0);
  PTR(lsst::daf::base::PropertySet) propSet = afwImg::readMetadata(fileName);
  PTR(afwImg::Wcs) wcs = afwImg::makeWcs(propSet);
  const afwImg::TanWcs* tanWcs = dynamic_cast<const afwImg::TanWcs*>(wcs.get());
  simAstrom::TanSipPix2RaDec gtransfoWcs = ConvertTanWcs(tanWcs);
  simAstrom::Point where(100.,200.);
  std::cout << std::setprecision(12) << "Poloka: " << gtransfoWcs.apply(where)<< std::endl;
  lsst::afw::geom::Point2D whereSame(100.,200.);
  PTR(lsst::afw::coord::Coord) coord = wcs->pixelToSky(whereSame);
  lsst::afw::geom::Point2D outDeg = coord->getPosition(lsst::afw::geom::degrees);
  std::cout << "Stack : " << outDeg[0] << ' ' << outDeg[1] << std::endl;
  return 0;

}
