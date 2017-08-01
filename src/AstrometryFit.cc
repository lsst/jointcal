#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "Eigen/Sparse"

#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/AstrometryFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Tripletlist.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.AstrometryFit");
}

namespace lsst {
namespace jointcal {

AstrometryFit::AstrometryFit(std::shared_ptr<Associations> associations,
                             std::shared_ptr<AstrometryModel> astrometryModel, double posError)
        : FitterBase(associations),
          _astrometryModel(astrometryModel),
          _refractionCoefficient(0),
          _nParDistortions(0),
          _nParPositions(0),
          _nParRefrac(_associations->getNFilters()),
          _posError(posError) {
    _JDRef = 0;

    _posError = posError;

    _referenceColor = 0;
    _sigCol = 0;
    unsigned count = 0;
    for (auto const &i : _associations->fittedStarList) {
        _referenceColor += i->color;
        _sigCol += std::pow(i->color, 2);
        count++;
    }
    if (count) {
        _referenceColor /= double(count);
        if (_sigCol > 0) _sigCol = sqrt(_sigCol / count - std::pow(_referenceColor, 2));
    }
    LOGLS_INFO(_log, "Reference Color: " << _referenceColor << " sig " << _sigCol);
}

#define NPAR_PM 2

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstrometryFit::transformFittedStar(FittedStar const &fittedStar, Gtransfo const *sky2TP,
                                         Point const &refractionVector, double refractionCoeff,
                                         double mjd) const {
    Point fittedStarInTP = sky2TP->apply(fittedStar);
    if (fittedStar.mightMove) {
        fittedStarInTP.x += fittedStar.pmx * mjd;
        fittedStarInTP.y += fittedStar.pmy * mjd;
    }
    // account for atmospheric refraction: does nothing if color
    // have not been assigned
    // the color definition shouldbe the same when computing derivatives
    double color = fittedStar.color - _referenceColor;
    fittedStarInTP.x += refractionVector.x * color * refractionCoeff;
    fittedStarInTP.y += refractionVector.y * color * refractionCoeff;
    return fittedStarInTP;
}

/*! This is the first implementation of an error "model".  We'll
  certainly have to upgrade it. MeasuredStar provides the mag in case
  we need it.  */
static void tweakAstromMeasurementErrors(FatPoint &P, MeasuredStar const &Ms, double error) {
    static bool called = false;
    static double increment = 0;
    if (!called) {
        increment = std::pow(error, 2);  // was in Preferences
        called = true;
    }
    P.vx += increment;
    P.vy += increment;
}

// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstrometryFit::leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                                      Eigen::VectorXd &fullGrad,
                                                      MeasuredStarList const *msList) const {
    /**********************************************************************/
    /* @note the math in this method and accumulateStatImage() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* Setup */
    /* this routine works in two different ways: either providing the
       ccdImage, of providing the MeasuredStarList. In the latter case, the
       ccdImage should match the one(s) in the list. */
    if (msList) assert(&(msList->front()->getCcdImage()) == &ccdImage);

    // get the Mapping
    const Mapping *mapping = _astrometryModel->getMapping(ccdImage);
    // count parameters
    unsigned npar_mapping = (_fittingDistortions) ? mapping->getNpar() : 0;
    unsigned npar_pos = (_fittingPos) ? 2 : 0;
    unsigned npar_refrac = (_fittingRefrac) ? 1 : 0;
    unsigned npar_pm = (_fittingPM) ? NPAR_PM : 0;
    unsigned npar_tot = npar_mapping + npar_pos + npar_refrac + npar_pm;
    // if (npar_tot == 0) this CcdImage does not contribute
    // any constraint to the fit, so :
    if (npar_tot == 0) return;
    std::vector<unsigned> indices(npar_tot, -1);
    if (_fittingDistortions) mapping->setMappingIndices(indices);

    // proper motion stuff
    double mjd = ccdImage.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = ccdImage.getRefractionVector();
    // transformation from sky to TP
    const Gtransfo *sky2TP = _astrometryModel->getSky2TP(ccdImage);
    // reserve matrices once for all measurements
    GtransfoLin dypdy;
    // the shape of H (et al) is required this way in order to be able to
    // separate derivatives along x and y as vectors.
    Eigen::MatrixX2d H(npar_tot, 2), halpha(npar_tot, 2), HW(npar_tot, 2);
    Eigen::Matrix2d transW(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::VectorXd grad(npar_tot);
    // current position in the Jacobian
    unsigned kTriplets = tripletList.getNextFreeIndex();
    const MeasuredStarList &catalog = (msList) ? *msList : ccdImage.getCatalogForFit();

    for (auto &i : catalog) {
        const MeasuredStar &ms = *i;
        if (!ms.isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = ms;
        tweakAstromMeasurementErrors(inPos, ms, _posError);
        H.setZero();  // we cannot be sure that all entries will be overwritten.
        FatPoint outPos;
        // should *not* fill H if whatToFit excludes mapping parameters.
        if (_fittingDistortions)
            mapping->computeTransformAndDerivatives(inPos, outPos, H);
        else
            mapping->transformPosAndErrors(inPos, outPos);

        unsigned ipar = npar_mapping;
        double det = outPos.vx * outPos.vy - std::pow(outPos.vxy, 2);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, "Inconsistent measurement errors: drop measurement at "
                                     << Point(ms) << " in image " << ccdImage.getName());
            continue;
        }
        transW(0, 0) = outPos.vy / det;
        transW(1, 1) = outPos.vx / det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy / det;
        // compute alpha, a triangular square root
        // of transW (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(transW(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = transW(0, 1) / alpha(0, 0);
        // DB - I think that the next line is equivalent to : alpha(1,1) = 1./sqrt(outPos.vy)
        // PA - seems correct !
        alpha(1, 1) = 1. / sqrt(det * transW(0, 0));
        alpha(0, 1) = 0;

        auto fs = ms.getFittedStar();

        Point fittedStarInTP =
                transformFittedStar(*fs, sky2TP, refractionVector, _refractionCoefficient, mjd);

        // compute derivative of TP position w.r.t sky position ....
        if (npar_pos > 0)  // ... if actually fitting FittedStar position
        {
            sky2TP->computeDerivative(*fs, dypdy, 1e-3);
            // sign checked
            // TODO Still have to check with non trivial non-diagonal terms
            H(npar_mapping, 0) = -dypdy.A11();
            H(npar_mapping + 1, 0) = -dypdy.A12();
            H(npar_mapping, 1) = -dypdy.A21();
            H(npar_mapping + 1, 1) = -dypdy.A22();
            indices[npar_mapping] = fs->getIndexInMatrix();
            indices.at(npar_mapping + 1) = fs->getIndexInMatrix() + 1;
            ipar += npar_pos;
        }
        /* only consider proper motions of objects allowed to move,
        unless the fit is going to be degenerate */
        if (_fittingPM && fs->mightMove) {
            H(ipar, 0) = -mjd;  // Sign unchecked but consistent with above
            H(ipar + 1, 1) = -mjd;
            indices[ipar] = fs->getIndexInMatrix() + 2;
            indices[ipar + 1] = fs->getIndexInMatrix() + 3;
            ipar += npar_pm;
        }
        if (_fittingRefrac) {
            /* if the definition of color changes, it has to remain
               consistent with transformFittedStar */
            double color = fs->color - _referenceColor;
            // sign checked
            H(ipar, 0) = -refractionVector.x * color;
            H(ipar, 1) = -refractionVector.y * color;
            indices[ipar] = _refracPosInMatrix;
            ipar += 1;
        }

        // We can now compute the residual
        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);

        // do not write grad = H*transW*res to avoid
        // dynamic allocation of a temporary
        halpha = H * alpha;
        HW = H * transW;
        grad = HW * res;
        // now feed in triplets and fullGrad
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar) {
            for (unsigned ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tripletList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            fullGrad(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }                    // end loop on measurements
    tripletList.setNextFreeIndex(kTriplets);
}

void AstrometryFit::leastSquareDerivativesReference(FittedStarList const &fittedStarList,
                                                    TripletList &tripletList,
                                                    Eigen::VectorXd &fullGrad) const {
    /**********************************************************************/
    /* @note the math in this method and accumulateStatRefStars() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* We compute here the derivatives of the terms involving fitted
       stars and reference stars. They only provide contributions if we
       are fitting positions: */
    if (!_fittingPos) return;
    /* the other case where the accumulation of derivatives stops
       here is when there are no RefStars */
    if (_associations->refStarList.size() == 0) return;
    Eigen::Matrix2d W(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::Matrix2d H(2, 2), halpha(2, 2), HW(2, 2);
    GtransfoLin der;
    Eigen::Vector2d res, grad;
    unsigned indices[2 + NPAR_PM];
    unsigned kTriplets = tripletList.getNextFreeIndex();
    /* We cannot use the spherical coordinates directly to evaluate
       Euclidean distances, we have to use a projector on some plane in
       order to express least squares. Not projecting could lead to a
       disaster around the poles or across alpha=0.  So we need a
       projector. We construct a projector and will change its
       projection point at every object */
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i : fittedStarList) {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        // Compute the derivative of the projector to incorporate its effects on the errors.
        proj.computeDerivative(fs, der, 1e-4);
        // sign checked. TODO check that the off-diagonal terms are OK.
        H(0, 0) = -der.A11();
        H(1, 0) = -der.A12();
        H(0, 1) = -der.A21();
        H(1, 1) = -der.A22();
        // TO DO : account for proper motions.
        double det = rsProj.vx * rsProj.vy - std::pow(rsProj.vxy, 2);
        if (rsProj.vx <= 0 || rsProj.vy <= 0 || det <= 0) {
            LOGLS_WARN(_log, "RefStar error matrix not positive definite for:  " << *rs);
            continue;
        }
        W(0, 0) = rsProj.vy / det;
        W(0, 1) = W(1, 0) = -rsProj.vxy / det;
        W(1, 1) = rsProj.vx / det;
        // compute alpha, a triangular square root
        // of W (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(W(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = W(0, 1) / alpha(0, 0);
        alpha(1, 1) = 1. / sqrt(det * W(0, 0));
        alpha(0, 1) = 0;
        indices[0] = fs.getIndexInMatrix();
        indices[1] = fs.getIndexInMatrix() + 1;
        unsigned npar_tot = 2;
        /* TODO: account here for proper motions in the reference
        catalog. We can code the effect and set the value to 0. Most
        (all?)  catalogs do not even come with a reference epoch. Gaia
        will change that. When refraction enters into the game, one should
        pay attention to the orientation of the frame */

        /* The residual should be Proj(fs)-Proj(*rs) in order to be consistent
        with the measurement terms. Since P(fs) = 0, we have: */
        res[0] = -rsProj.x;
        res[1] = -rsProj.y;
        halpha = H * alpha;
        // grad = H*W*res
        HW = H * W;
        grad = HW * res;
        // now feed in triplets and fullGrad
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar) {
            for (unsigned ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tripletList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            fullGrad(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }
    tripletList.setNextFreeIndex(kTriplets);
}

void AstrometryFit::accumulateStatImage(CcdImage const &ccdImage, Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesMeasurement() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/
    /* Setup */
    // 1 : get the Mapping's
    const Mapping *mapping = _astrometryModel->getMapping(ccdImage);
    // proper motion stuff
    double mjd = ccdImage.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = ccdImage.getRefractionVector();
    // transformation from sky to TP
    const Gtransfo *sky2TP = _astrometryModel->getSky2TP(ccdImage);
    // reserve matrix once for all measurements
    Eigen::Matrix2Xd transW(2, 2);

    auto &catalog = ccdImage.getCatalogForFit();
    for (auto const &ms : catalog) {
        if (!ms->isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = *ms;
        tweakAstromMeasurementErrors(inPos, *ms, _posError);

        FatPoint outPos;
        // should *not* fill H if whatToFit excludes mapping parameters.
        mapping->transformPosAndErrors(inPos, outPos);
        double det = outPos.vx * outPos.vy - std::pow(outPos.vxy, 2);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, " Inconsistent measurement errors :drop measurement at "
                                     << Point(*ms) << " in image " << ccdImage.getName());
            continue;
        }
        transW(0, 0) = outPos.vy / det;
        transW(1, 1) = outPos.vx / det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy / det;

        auto fs = ms->getFittedStar();
        Point fittedStarInTP =
                transformFittedStar(*fs, sky2TP, refractionVector, _refractionCoefficient, mjd);

        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);
        double chi2Val = res.transpose() * transW * res;

        accum.addEntry(chi2Val, 2, ms);
    }  // end of loop on measurements
}

void AstrometryFit::accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const {
    for (auto const &ccdImage : ccdImageList) {
        accumulateStatImage(*ccdImage, accum);
    }
}

void AstrometryFit::accumulateStatRefStars(Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesReference() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* If you wonder why we project here, read comments in
       AstrometryFit::leastSquareDerivativesReference(TripletList &TList, Eigen::VectorXd &Rhs) */
    FittedStarList &fittedStarList = _associations->fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &fs : fittedStarList) {
        const RefStar *rs = fs->getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(*fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        // TO DO : account for proper motions.
        double rx = rsProj.x;  // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx * rsProj.vy - std::pow(rsProj.vxy, 2);
        double wxx = rsProj.vy / det;
        double wyy = rsProj.vx / det;
        double wxy = -rsProj.vxy / det;
        accum.addEntry(wxx * std::pow(rx, 2) + 2 * wxy * rx * ry + wyy * std::pow(ry, 2), 2, fs);
    }
}

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void AstrometryFit::getIndicesOfMeasuredStar(MeasuredStar const &measuredStar,
                                             std::vector<unsigned> &indices) const {
    if (_fittingDistortions) {
        const Mapping *mapping = _astrometryModel->getMapping(measuredStar.getCcdImage());
        mapping->setMappingIndices(indices);
    }
    auto fs = measuredStar.getFittedStar();
    unsigned fsIndex = fs->getIndexInMatrix();
    if (_fittingPos) {
        indices.push_back(fsIndex);
        indices.push_back(fsIndex + 1);
    }
    // For securing the outlier removal, the next block is just useless
    if (_fittingPM) {
        for (unsigned k = 0; k < NPAR_PM; ++k) indices.push_back(fsIndex + 2 + k);
    }
    /* Should not put the index of refaction stuff or we will not be
       able to remove more than 1 star at a time. */
}

void AstrometryFit::assignIndices(std::string const &whatToFit) {
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: Now fitting " << whatToFit);
    _fittingDistortions = (_whatToFit.find("Distortions") != std::string::npos);
    _fittingPos = (_whatToFit.find("Positions") != std::string::npos);
    _fittingRefrac = (_whatToFit.find("Refrac") != std::string::npos);
    if (_sigCol == 0 && _fittingRefrac) {
        LOGLS_WARN(_log,
                   "Cannot fit refraction coefficients without a color lever arm. Ignoring refraction.");
        _fittingRefrac = false;
    }
    _fittingPM = (_whatToFit.find("PM") != std::string::npos);
    // When entering here, we assume that whatToFit has already been interpreted.

    _nParDistortions = 0;
    if (_fittingDistortions) _nParDistortions = _astrometryModel->assignIndices(0, _whatToFit);
    unsigned ipar = _nParDistortions;

    if (_fittingPos) {
        FittedStarList &fittedStarList = _associations->fittedStarList;
        for (auto &fittedStar : fittedStarList) {
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in GetMeasuredStarIndices
            fittedStar->setIndexInMatrix(ipar);
            ipar += 2;
            if ((_fittingPM)&fittedStar->mightMove) ipar += NPAR_PM;
        }
    }
    _nParPositions = ipar - _nParDistortions;
    if (_fittingRefrac) {
        _refracPosInMatrix = ipar;
        ipar += _nParRefrac;
    }
    _nParTot = ipar;
}

void AstrometryFit::offsetParams(Eigen::VectorXd const &delta) {
    if (delta.size() != _nParTot)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "AstrometryFit::offsetParams : the provided vector length is not compatible with "
                          "the current whatToFit setting");
    if (_fittingDistortions) _astrometryModel->offsetParams(delta);

    if (_fittingPos) {
        FittedStarList &fittedStarList = _associations->fittedStarList;
        for (auto const &i : fittedStarList) {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when assigning indices (assignIndices())
            unsigned index = fs.getIndexInMatrix();
            fs.x += delta(index);
            fs.y += delta(index + 1);
            if ((_fittingPM)&fs.mightMove) {
                fs.pmx += delta(index + 2);
                fs.pmy += delta(index + 3);
            }
        }
    }
    if (_fittingRefrac) {
        _refractionCoefficient += delta(_refracPosInMatrix);
    }
}

// should not be too large !
#ifdef STORAGE
static void write_sparse_matrix_in_fits(SpMat const &mat, std::string const &fitsName) {
    if (mat.rows() * mat.cols() > 2e8) {
        LOGLS_WARN(_log,
                   "write_sparse_matrix_in_fits: yout matrix is too large. " << fitsName << " not generated");
        return;
    }
    Mat m(mat.rows(), mat.cols());
    for (int k = 0; k < mat.outerSize(); ++k)
        for (SpMat::InnerIterator it(mat, k); it; ++it) {
            m(it.row(), it.col()) = it.value();
        }
    m.writeFits(fitsName);
}

static void write_vect_in_fits(Eigen::VectorXd const &vectorXd, std::string const &fitsName) {
    Vect v(vectorXd.size());
    for (int k = 0; k < vectorXd.size(); ++k) v(k) = V(k);
    Mat(v).writeFits(fitsName);
}

#endif

void AstrometryFit::checkStuff() {
#if (0)
    const char *what2fit[] = {"Positions",
                              "Distortions",
                              "Refrac",
                              "Positions Distortions",
                              "Positions Refrac",
                              "Distortions Refrac",
                              "Positions Distortions Refrac"};
#endif
    const char *what2fit[] = {"Positions", "Distortions", "Positions Distortions"};
    // DEBUG
    for (unsigned k = 0; k < sizeof(what2fit) / sizeof(what2fit[0]); ++k) {
        assignIndices(what2fit[k]);
        TripletList tripletList(10000);
        Eigen::VectorXd grad(_nParTot);
        grad.setZero();
        leastSquareDerivatives(tripletList, grad);
        SpMat jacobian(_nParTot, tripletList.getNextFreeIndex());
        jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
        SpMat hessian = jacobian * jacobian.transpose();
#ifdef STORAGE
        char name[24];
        sprintf(name, "H%d.fits", k);
        write_sparse_matrix_in_fits(hessian, name);
        sprintf(name, "g%d.fits", k);
        write_vect_in_fits(grad, name);
#endif
        LOGLS_DEBUG(_log, "npar : " << _nParTot << ' ' << _nParDistortions);
    }
}

void AstrometryFit::saveResultTuples(std::string const &tupleName) const {
    /* cook-up 2 different file names by inserting something just before
       the dot (if any), and within the actual file name. */
    size_t dot = tupleName.rfind('.');
    size_t slash = tupleName.rfind('/');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) dot = tupleName.size();
    std::string meas_tuple(tupleName);
    meas_tuple.insert(dot, "-meas");
    makeMeasResTuple(meas_tuple);
    std::string ref_tuple(tupleName);
    ref_tuple.insert(dot, "-ref");
    makeRefResTuple(ref_tuple);
}

void AstrometryFit::makeMeasResTuple(std::string const &tupleName) const {
    std::ofstream tuple(tupleName.c_str());
    tuple << "#xccd: coordinate in CCD" << std::endl
          << "#yccd: " << std::endl
          << "#rx:   residual in degrees in TP" << std::endl
          << "#ry:" << std::endl
          << "#xtp: transformed coordinate in TP " << std::endl
          << "#ytp:" << std::endl
          << "#mag: rough mag" << std::endl
          << "#jd: Julian date of the measurement" << std::endl
          << "#rvx: transformed measurement uncertainty " << std::endl
          << "#rvy:" << std::endl
          << "#rvxy:" << std::endl
          << "#color : " << std::endl
          << "#fsindex: some unique index of the object" << std::endl
          << "#ra: pos of fitted star" << std::endl
          << "#dec: pos of fitted star" << std::endl
          << "#chi2: contribution to Chi2 (2D dofs)" << std::endl
          << "#nm: number of measurements of this FittedStar" << std::endl
          << "#chip: chip number" << std::endl
          << "#visit: visit id" << std::endl
          << "#end" << std::endl;
    const CcdImageList &L = _associations->getCcdImageList();
    for (auto const &i : L) {
        const CcdImage &im = *i;
        const MeasuredStarList &cat = im.getCatalogForFit();
        const Mapping *mapping = _astrometryModel->getMapping(im);
        const Point &refractionVector = im.getRefractionVector();
        double mjd = im.getMjd() - _JDRef;
        for (auto const &is : cat) {
            const MeasuredStar &ms = *is;
            if (!ms.isValid()) continue;
            FatPoint tpPos;
            FatPoint inPos = ms;
            tweakAstromMeasurementErrors(inPos, ms, _posError);
            mapping->transformPosAndErrors(inPos, tpPos);
            const Gtransfo *sky2TP = _astrometryModel->getSky2TP(im);
            auto fs = ms.getFittedStar();

            Point fittedStarInTP =
                    transformFittedStar(*fs, sky2TP, refractionVector, _refractionCoefficient, mjd);
            Point res = tpPos - fittedStarInTP;
            double det = tpPos.vx * tpPos.vy - std::pow(tpPos.vxy, 2);
            double wxx = tpPos.vy / det;
            double wyy = tpPos.vx / det;
            double wxy = -tpPos.vxy / det;
            //      double chi2 = rx*(wxx*rx+wxy*ry)+ry*(wxy*rx+wyy*ry);
            double chi2 = wxx * res.x * res.x + wyy * res.y * res.y + 2 * wxy * res.x * res.y;
            tuple << std::setprecision(9);
            tuple << ms.x << ' ' << ms.y << ' ' << res.x << ' ' << res.y << ' ' << tpPos.x << ' ' << tpPos.y
                  << ' ' << fs->getMag() << ' ' << mjd << ' ' << tpPos.vx << ' ' << tpPos.vy << ' '
                  << tpPos.vxy << ' ' << fs->color << ' ' << fs->getIndexInMatrix() << ' ' << fs->x << ' '
                  << fs->y << ' ' << chi2 << ' ' << fs->getMeasurementCount() << ' ' << im.getCcdId() << ' '
                  << im.getVisit() << std::endl;
        }  // loop on measurements in image
    }      // loop on images
}

void AstrometryFit::makeRefResTuple(std::string const &tupleName) const {
    std::ofstream tuple(tupleName.c_str());
    tuple << "#ra: coordinates of FittedStar" << std::endl
          << "#dec: " << std::endl
          << "#rx:   residual in degrees in TP" << std::endl
          << "#ry:" << std::endl
          << "#mag: mag" << std::endl
          << "#rvx: transformed measurement uncertainty " << std::endl
          << "#rvy:" << std::endl
          << "#rvxy:" << std::endl
          << "#color : " << std::endl
          << "#fsindex: some unique index of the object" << std::endl
          << "#chi2: contribution to Chi2 (2D dofs)" << std::endl
          << "#nm: number of measurements of this FittedStar" << std::endl
          << "#end" << std::endl;
    // The following loop is heavily inspired from AstrometryFit::computeChi2()
    const FittedStarList &fittedStarList = _associations->fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i : fittedStarList) {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        double rx = rsProj.x;  // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx * rsProj.vy - std::pow(rsProj.vxy, 2);
        double wxx = rsProj.vy / det;
        double wyy = rsProj.vx / det;
        double wxy = -rsProj.vxy / det;
        double chi2 = wxx * std::pow(rx, 2) + 2 * wxy * rx * ry + wyy * std::pow(ry, 2);
        tuple << std::setprecision(9);
        tuple << fs.x << ' ' << fs.y << ' ' << rx << ' ' << ry << ' ' << fs.getMag() << ' ' << rsProj.vx
              << ' ' << rsProj.vy << ' ' << rsProj.vxy << ' ' << fs.color << ' ' << fs.getIndexInMatrix()
              << ' ' << chi2 << ' ' << fs.getMeasurementCount() << std::endl;
    }  // loop on FittedStars
}
}  // namespace jointcal
}  // namespace lsst
