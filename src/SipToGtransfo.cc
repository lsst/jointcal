#include <algorithm>

#include "Eigen/Core"
#include "lsst/meas/simastrom/SipToGtransfo.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/meas/simastrom/Point.h"
#include "lsst/meas/simastrom/Frame.h"
#include "lsst/daf/base/PropertySet.h"

namespace simAstrom = lsst::meas::simastrom; 
namespace afwImg = lsst::afw::image;
namespace afwGeom  = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace simastrom {

static const int lsstToFitsPixels = +1;
static const int fitsToLsstPixels = -1;

typedef boost::shared_ptr<simAstrom::GtransfoPoly> GtPoly_Ptr;
	

simAstrom::TanSipPix2RaDec ConvertTanWcs(const boost::shared_ptr<lsst::afw::image::TanWcs> wcs)
{
  GtPoly_Ptr sipCorr(new simAstrom::GtransfoPoly(0));

  /* TanWcs::pixelToSkyImpl adds  - PixelZeroPos + lsstToFitsPixels
     to the input coordinates. */

  simAstrom::GtransfoLinShift firstShift(-afwImg::PixelZeroPos + lsstToFitsPixels,					 -afwImg::PixelZeroPos + lsstToFitsPixels);

  /* beware : Wcs::getPixelOrigin return crpix + fitsToLsstPixels */
  lsst::afw::geom::Point2D offset_crpix = wcs->getPixelOrigin();

  lsst::daf::base::PropertyList::Ptr wcsMeta = wcs->getFitsMetadata();

  if (wcs->hasDistortion())
    {
      Eigen::MatrixXd sipA;
      Eigen::MatrixXd sipB;
      lsst::afw::image::TanWcs::decodeSipHeader(*wcsMeta, "A", sipA);
      lsst::afw::image::TanWcs::decodeSipHeader(*wcsMeta, "B", sipB);
  
      int sipOrder = std::max(wcsMeta->get<int>("A_ORDER"), wcsMeta->get<int>("B_ORDER"));

      simAstrom::GtransfoPoly sipPoly(sipOrder);
      for (int i=0; i<=sipOrder; ++i)
        {
	  for (int j=0; j<=sipOrder; ++j)
	    {
	      if (i<sipA.cols() && j<sipA.rows()&& (i+j)<= sipOrder) 
	        sipPoly.Coeff(i,j,0) = sipA(i,j);
	      if (i<sipB.cols() && j<sipB.rows() && (i+j)<= sipOrder) 
	        sipPoly.Coeff(i,j,1) = sipB(i,j);
	    }
	}

      /*
         TanWcs::undistortPixel adds (- wcs_info->crpix) (line 330) to
	 the coordinates before computing the polynomials so the total
	 shift reads:

	 -PixelZeroPos + lsstToFitsPixels -(offset_crpix-fitsToLsstPixels)
	 = -PixelZeroPos -offset_crpix

	 remember that offset_crpix = crpix + fitsToLsstPixels,
	 so crpix = offset_crpix - fitsToLsstPixels,
	 which is the  second coordinate shift:
      */

      simAstrom::GtransfoLinShift secondShift(-offset_crpix[0]+fitsToLsstPixels,
					      -offset_crpix[1]+fitsToLsstPixels);

      /* then the SIP correction (TanWcs::undistorPixel, last line)
	 returns pix + sipPoly*secondShift(pix) where pix =
	 firstShift(input) */

      simAstrom::GtransfoPoly actualSip = firstShift+
	sipPoly*secondShift*firstShift;
      	
      sipCorr.reset(new simAstrom::GtransfoPoly(actualSip));
    }

  // now compute the lin part (nothing to do with SIP) */
  Eigen::Matrix2d cdMat = wcs->getCDMatrix();
  simAstrom::GtransfoLin cdTrans;
  cdTrans.Coeff(1,0,0) = cdMat(0,0); // CD1_1
  cdTrans.Coeff(0,1,0) = cdMat(0,1); // CD1_2 
  cdTrans.Coeff(1,0,1) = cdMat(1,0); // CD2_1
  cdTrans.Coeff(0,1,1) = cdMat(1,1); // CD2_1
  // this is by chance equal to secondShift. We do not rely on that.
  simAstrom::GtransfoLinShift crpixShift(-offset_crpix[0] +  fitsToLsstPixels,
					 -offset_crpix[1] +  fitsToLsstPixels); 

  // CD's apply to CRPIX-shifted coordinate
  simAstrom::GtransfoLin linPart = cdTrans * crpixShift;

  /* It there are no distorsions, we have to apply firstShift in order
     to emulateTanWcs::pixelToSkyImpl. If there are distorsions, this step
     is already included in the sip corrections. */
  if (!wcs->hasDistortion()) 
    linPart = linPart* firstShift;

  //  lsst::afw::coord::Coord tp = wcs->getSkyOrigin()->getPosition(lsst::afw::geom::degrees);
  // the above line returns radians ?!
  double ra  = wcsMeta->get<double>("CRVAL1");
  double dec = wcsMeta->get<double>("CRVAL2");
  
  simAstrom::Point tangentPoint(ra,dec);

  // return simAstrom::TanSipPix2RaDec(linPart, tangentPoint, sipCorr->get());
  return simAstrom::TanSipPix2RaDec(linPart, tangentPoint, (const simAstrom::GtransfoPoly*) sipCorr.get());

}



/* The inverse transformation i.e. convert from the fit result to the SIP
   convention. */
PTR(afwImg::TanWcs) GtransfoToSip(const simAstrom::TanSipPix2RaDec WcsTransfo, 
				  const simAstrom::Frame &CcdFrame)
{
  GtransfoLin linPart = WcsTransfo.LinPart();
  afwGeom::Point2D crpix_lsst; // in LSST "frame"
  // compute crpix as the point that transforms to (0,0):
  linPart.invert().apply(0.,0., crpix_lsst[0], crpix_lsst[1]);
  /* At this stage, crpix should not be shifted from "LSST units" to
     "FITS units" yet because the TanWcs constructors expect it in LSST
     units */


  // crval from type conversion
  afwGeom::Point2D crval;
  crval[0] = WcsTransfo.TangentPoint().x;
  crval[1] = WcsTransfo.TangentPoint().y;

  // CD matrix:
  Eigen::Matrix2d cdMat;
  cdMat(0,0) = linPart.Coeff(1,0,0); // CD1_1
  cdMat(0,1) = linPart.Coeff(0,1,0); // CD1_2 
  cdMat(1,0) = linPart.Coeff(1,0,1); // CD2_1
  cdMat(1,1) = linPart.Coeff(0,1,1); // CD2_1
  
  if (!WcsTransfo.Corr()) // the WCS has no distortions
    return boost::shared_ptr<afwImg::TanWcs>(new afwImg::TanWcs(crval,crpix_lsst,cdMat));

  /* The CRPIX in the fits header (and hence applied by the wcslib)
     will be different from the one we computed above */
  double offset = -afwImg::PixelZeroPos+lsstToFitsPixels;
  simAstrom::GtransfoLinShift s1(offset,offset);
  linPart = linPart*s1.invert(); // this is what the wcslib will apply
  
  // This is (the opposite of) the crpix that will go into the fits header:
  simAstrom::GtransfoLinShift s2(-crpix_lsst[0]-offset,-crpix_lsst[1]-offset);


  // for SIP, pix2TP = linpart*sip, so
  simAstrom::GtransfoPoly sipTransform = simAstrom::GtransfoPoly(linPart.invert())*WcsTransfo.Pix2TangentPlane();
  // then the sip transform reads ST = (ID+PA*S2)*S1, so PA =(ST*S1^-1-ID)*S2^-1 
  simAstrom::GtransfoLin id; // default constructor
  simAstrom::GtransfoPoly sipPoly = (sipTransform*s1.invert()-id)*s2.invert();

  // there is still a 1 pixel offset somewhere. To be found.

  int sipOrder = sipPoly.Degree();
  Eigen::MatrixXd sipA(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  Eigen::MatrixXd sipB(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  for (int i=0; i<=sipOrder; ++i)
    for (int j=0; j<=sipOrder-i; ++j)
      {
	sipA(i,j) = sipPoly.Coeff(i,j,0);
	sipB(i,j) = sipPoly.Coeff(i,j,1);
      }

  // No reverse stuff for now
  Eigen::MatrixXd sipAp(Eigen::MatrixXd::Zero(sipOrder,sipOrder));
  Eigen::MatrixXd sipBp(Eigen::MatrixXd::Zero(sipOrder,sipOrder));

  return boost::shared_ptr<afwImg::TanWcs>(new afwImg::TanWcs(crval, crpix_lsst, cdMat, sipA, sipB, sipAp, sipBp));
  
  /* now, one subtelty of SIP: the direct (A,B) and inverse(Ap, Bp) 
     polynomials are not inverse of each other */
  
  

  






}

}}} 

  
   



  

