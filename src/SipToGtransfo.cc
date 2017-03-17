#include <algorithm>

#include "Eigen/Core"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/daf/base/PropertySet.h"

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;
namespace afwGeom  = lsst::afw::geom;

namespace lsst {
namespace jointcal {

static const int lsstToFitsPixels = +1;
static const int fitsToLsstPixels = -1;

typedef std::shared_ptr<jointcal::GtransfoPoly> GtPoly_Ptr;


jointcal::TanSipPix2RaDec convertTanWcs(const std::shared_ptr<lsst::afw::image::TanWcs> wcs)
{
  GtPoly_Ptr sipCorr(new jointcal::GtransfoPoly(0));

  /* beware : Wcs::getPixelOrigin return crpix_fits +
   fitsToLsstPixels, so all the algebra we perform here happens in the
   "Lsst frame", i.e (0,0)-based. this algebra is justified in the
   documentation of the package. */

  lsst::afw::geom::Point2D crpix_lsst = wcs->getPixelOrigin();

  lsst::daf::base::PropertyList::Ptr wcsMeta = wcs->getFitsMetadata();

  if (wcs->hasDistortion())
    {
      Eigen::MatrixXd sipA;
      Eigen::MatrixXd sipB;
      lsst::afw::image::TanWcs::decodeSipHeader(*wcsMeta, "A", sipA);
      lsst::afw::image::TanWcs::decodeSipHeader(*wcsMeta, "B", sipB);

      int sipOrder = std::max(wcsMeta->get<int>("A_ORDER"), wcsMeta->get<int>("B_ORDER"));

      jointcal::GtransfoPoly sipPoly(sipOrder);
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

      jointcal::GtransfoLinShift s2(-crpix_lsst[0], - crpix_lsst[1]);

      /* then the SIP correction (TanWcs::undistorPixel, last line)
	 returns pix + sipPoly*secondShift(pix) where secondShift
	 subtracts crpix_header from (1,1) based coordinates, i.e. the
	 same thing as subtracting crpix_lsst from (0,0-based
	 coordinates. So undistort pixel does:
            id+sipPoly*s2
      */

      GtransfoLin id; // identity is the default constructor.
      // This is what is returned by TanWcs::undistortpixel
      jointcal::GtransfoPoly actualSip = id + sipPoly*s2;

      sipCorr.reset(new jointcal::GtransfoPoly(actualSip));
    }

  // now compute the lin part (nothing to do with SIP) */
  Eigen::Matrix2d cdMat = wcs->getCDMatrix();
  jointcal::GtransfoLin cdTrans;
  cdTrans.Coeff(1,0,0) = cdMat(0,0); // CD1_1
  cdTrans.Coeff(0,1,0) = cdMat(0,1); // CD1_2
  cdTrans.Coeff(1,0,1) = cdMat(1,0); // CD2_1
  cdTrans.Coeff(0,1,1) = cdMat(1,1); // CD2_1
  // this is by chance equal to s2, but we will not rely on this fact:
  jointcal::GtransfoLinShift crpixShift(-crpix_lsst[0], -crpix_lsst[1]);

  // CD's apply to CRPIX-shifted coordinate
  jointcal::GtransfoLin linPart = cdTrans * crpixShift;

  //  lsst::afw::coord::Coord tp = wcs->getSkyOrigin()->getPosition(lsst::afw::geom::degrees);
  // the above line returns radians ?!
  double ra  = wcsMeta->get<double>("CRVAL1");
  double dec = wcsMeta->get<double>("CRVAL2");

  jointcal::Point tangentPoint(ra,dec);

  // return jointcal::TanSipPix2RaDec(linPart, tangentPoint, sipCorr->get());
  return jointcal::TanSipPix2RaDec(linPart, tangentPoint, (const jointcal::GtransfoPoly*) sipCorr.get());

}



/* The inverse transformation i.e. convert from the fit result to the SIP
   convention. */
PTR(afwImg::TanWcs) gtransfoToTanWcs(const jointcal::TanSipPix2RaDec WcsTransfo,
				     const jointcal::Frame &CcdFrame,
				     const bool NoLowOrderSipTerms)
{
  GtransfoLin linPart = WcsTransfo.LinPart();
  afwGeom::Point2D crpix_lsst; // in LSST "frame"
  /* In order to remove the low order sip terms, one has to
     define the linear WCS transformation as the expansion of
     the total pix-to-tangent plane (or focal plane) at the
     tangent point. In order to do that, we first have to find
     which pixel transforms to the tangent point, and then expand there */

  /* compute crpix as the point that is transformed by the linear part
    into (0,0) */
  linPart.invert().apply(0.,0., crpix_lsst[0], crpix_lsst[1]);

  // This is what we have to respect:
  jointcal::GtransfoPoly pix2TP = WcsTransfo.Pix2TangentPlane();

  if (NoLowOrderSipTerms) {
    Point ctmp = Point(crpix_lsst[0], crpix_lsst[1]);
    // cookup a large Frame
    jointcal::Frame f(ctmp.x-10000, ctmp.y-10000, ctmp.x+10000, ctmp.y +10000);
    jointcal::Gtransfo *r =pix2TP.InverseTransfo(1e-6, f);
    // overwrite crpix ...
    r->apply(0,0, crpix_lsst[0], crpix_lsst[1]);
    delete r;
    // and the "linpart"
    linPart = pix2TP.LinearApproximation(Point(crpix_lsst[0], crpix_lsst[1]));
  }

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
  cdMat(1,1) = linPart.Coeff(0,1,1); // CD2_2

  if (!WcsTransfo.Corr()) // the WCS has no distortions
    return std::shared_ptr<afwImg::TanWcs>(new afwImg::TanWcs(crval,crpix_lsst,cdMat));

  /* We are now given:
     - CRPIX
     - The CD matrix
     - the lin part (i.e. the combination of CRPIX and CD)
     - and pix2TP, the total transformation from pixels to tangent plane
     and we want to extract the SIP polynomials. The algebra is detailed
     in the appendix of the documentation */

  // This is (the opposite of) the crpix that will go into the fits header:
  jointcal::GtransfoLinShift s2(-crpix_lsst[0],-crpix_lsst[1]);

  // for SIP, pix2TP = linpart*sipStuff, so
  jointcal::GtransfoPoly sipTransform = jointcal::GtransfoPoly(linPart.invert())* pix2TP;
  // then the sip transform reads ST = (ID+PA*S2)
  //   PA*S2 = ST -ID,   PA = (ST-Id)*S2^-1
  jointcal::GtransfoLin id; // default constructor = identity
  jointcal::GtransfoPoly sipPoly = (sipTransform-id)*s2.invert();

  // coockup the inverse sip polynomials
  // last argument : precision in pixels.
  jointcal::GtransfoPoly *tp2Pix = InversePolyTransfo(pix2TP, CcdFrame, 1e-4);
  if (!tp2Pix)
    {
      throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "GtransfoToSip: could not invert the input wcs ");
    }
  jointcal::GtransfoPoly invSipStuff = (*tp2Pix)*linPart;
  delete tp2Pix;
  jointcal::GtransfoPoly sipPolyInv =  (invSipStuff -id)*s2.invert();

  // now extract sip coefficients. First forward ones:
  int sipOrder = sipPoly.Degree();
  Eigen::MatrixXd sipA(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  Eigen::MatrixXd sipB(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  for (int i=0; i<=sipOrder; ++i)
    for (int j=0; j<=sipOrder-i; ++j)
      {
	sipA(i,j) = sipPoly.Coeff(i,j,0);
	sipB(i,j) = sipPoly.Coeff(i,j,1);
      }

  // now backwards coefficients
  sipOrder = sipPolyInv.Degree();
  Eigen::MatrixXd sipAp(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  Eigen::MatrixXd sipBp(Eigen::MatrixXd::Zero(sipOrder+1,sipOrder+1));
  for (int i=0; i<=sipOrder; ++i)
    for (int j=0; j<=sipOrder-i; ++j)
      {
	sipAp(i,j) = sipPolyInv.Coeff(i,j,0);
	sipBp(i,j) = sipPolyInv.Coeff(i,j,1);
      }

  return std::shared_ptr<afwImg::TanWcs>(new afwImg::TanWcs(crval, crpix_lsst, cdMat, sipA, sipB, sipAp, sipBp));

}

}} // end of namespaces








