#include <algorithm>

#include "Eigen/Core"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/geom/SpherePoint.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/daf/base/PropertySet.h"

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace jointcal {

typedef std::shared_ptr<jointcal::GtransfoPoly> GtPoly_Ptr;

/* The inverse transformation i.e. convert from the fit result to the SIP
   convention. */
std::shared_ptr<afw::geom::SkyWcs> gtransfoToTanWcs(const jointcal::TanSipPix2RaDec wcsTransfo,
                                                    const jointcal::Frame &ccdFrame,
                                                    const bool noLowOrderSipTerms) {
    GtransfoLin linPart = wcsTransfo.getLinPart();
    afwGeom::Point2D crpixLsst;  // in LSST "frame"
                                  /* In order to remove the low order sip terms, one has to
                                     define the linear WCS transformation as the expansion of
                                     the total pix-to-tangent plane (or focal plane) at the
                                     tangent point. In order to do that, we first have to find
                                     which pixel transforms to the tangent point, and then expand there */

    /* compute crpix as the point that is transformed by the linear part
      into (0,0) */
    linPart.invert().apply(0., 0., crpixLsst[0], crpixLsst[1]);

    // This is what we have to respect:
    jointcal::GtransfoPoly pix2TP = wcsTransfo.getPix2TangentPlane();

    if (noLowOrderSipTerms) {
        Point ctmp = Point(crpixLsst[0], crpixLsst[1]);
        // cookup a large Frame
        jointcal::Frame f(ctmp.x - 10000, ctmp.y - 10000, ctmp.x + 10000, ctmp.y + 10000);
        auto r = pix2TP.inverseTransfo(1e-6, f);
        // overwrite crpix ...
        r->apply(0, 0, crpixLsst[0], crpixLsst[1]);
        // and the "linpart"
        linPart = pix2TP.linearApproximation(Point(crpixLsst[0], crpixLsst[1]));
    }

    /* At this stage, crpix should not be shifted from "LSST units" to
       "FITS units" yet because the SkyWcs constructors expect it in LSST
       units */

    afw::geom::SpherePoint const crval(wcsTransfo.getTangentPoint().x * afwGeom::DEGREES,
                                       wcsTransfo.getTangentPoint().y * afwGeom::DEGREES);

    // CD matrix:
    Eigen::Matrix2d cdMat;
    cdMat(0, 0) = linPart.coeff(1, 0, 0);  // CD1_1
    cdMat(0, 1) = linPart.coeff(0, 1, 0);  // CD1_2
    cdMat(1, 0) = linPart.coeff(1, 0, 1);  // CD2_1
    cdMat(1, 1) = linPart.coeff(0, 1, 1);  // CD2_2

    if (!wcsTransfo.getCorr()) {  // the WCS has no distortions
        return afw::geom::makeSkyWcs(crpixLsst, crval, cdMat);
    }

    /* We are now given:
       - CRPIX
       - The CD matrix
       - the lin part (i.e. the combination of CRPIX and CD)
       - and pix2TP, the total transformation from pixels to tangent plane
       and we want to extract the SIP polynomials. The algebra is detailed
       in the appendix of the documentation */

    // This is (the opposite of) the crpix that will go into the fits header:
    jointcal::GtransfoLinShift s2(-crpixLsst[0], -crpixLsst[1]);

    // for SIP, pix2TP = linpart*sipStuff, so
    jointcal::GtransfoPoly sipTransform = jointcal::GtransfoPoly(linPart.invert()) * pix2TP;
    // then the sip transform reads ST = (ID+PA*S2)
    //   PA*S2 = ST -ID,   PA = (ST-Id)*S2^-1
    jointcal::GtransfoLin id;  // default constructor = identity
    jointcal::GtransfoPoly sipPoly = (sipTransform - id) * s2.invert();

    // coockup the inverse sip polynomials
    // last argument : precision in pixels.
    auto tp2Pix = inversePolyTransfo(pix2TP, ccdFrame, 1e-4);
    if (!tp2Pix) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "GtransfoToSip: could not invert the input wcs ");
    }
    jointcal::GtransfoPoly invSipStuff = (*tp2Pix) * linPart;
    jointcal::GtransfoPoly sipPolyInv = (invSipStuff - id) * s2.invert();

    // now extract sip coefficients. First forward ones:
    int sipOrder = sipPoly.getDegree();
    Eigen::MatrixXd sipA(Eigen::MatrixXd::Zero(sipOrder + 1, sipOrder + 1));
    Eigen::MatrixXd sipB(Eigen::MatrixXd::Zero(sipOrder + 1, sipOrder + 1));
    for (int i = 0; i <= sipOrder; ++i) {
        for (int j = 0; j <= sipOrder - i; ++j) {
            sipA(i, j) = sipPoly.coeff(i, j, 0);
            sipB(i, j) = sipPoly.coeff(i, j, 1);
        }
    }

    // now backwards coefficients
    sipOrder = sipPolyInv.getDegree();
    Eigen::MatrixXd sipAp(Eigen::MatrixXd::Zero(sipOrder + 1, sipOrder + 1));
    Eigen::MatrixXd sipBp(Eigen::MatrixXd::Zero(sipOrder + 1, sipOrder + 1));
    for (int i = 0; i <= sipOrder; ++i) {
        for (int j = 0; j <= sipOrder - i; ++j) {
            sipAp(i, j) = sipPolyInv.coeff(i, j, 0);
            sipBp(i, j) = sipPolyInv.coeff(i, j, 1);
        }
    }

    return makeTanSipWcs(crpixLsst, crval, cdMat, sipA, sipB, sipAp, sipBp);
}

}  // namespace jointcal
}  // namespace lsst
