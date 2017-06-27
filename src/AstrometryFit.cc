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
#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Tripletlist.h"

using namespace std;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.AstrometryFit");
}

namespace lsst {
namespace jointcal {

AstrometryFit::AstrometryFit(Associations &associations, AstrometryModel *astrometryModel, double posError)
        : _associations(associations), _astrometryModel(astrometryModel), _posError(posError) {
    _LastNTrip = 0;
    _JDRef = 0;

    _posError = posError;

    _referenceColor = 0;
    _sigCol = 0;
    unsigned count = 0;
    for (auto const &i : _associations.fittedStarList) {
        _referenceColor += i->color;
        _sigCol += sqr(i->color);
        count++;
    }
    if (count) {
        _referenceColor /= double(count);
        if (_sigCol > 0) _sigCol = sqrt(_sigCol / count - sqr(_referenceColor));
    }
    LOGLS_INFO(_log, "Reference Color: " << _referenceColor << " sig " << _sigCol);

    _nRefrac = _associations.getNFilters();
    _refractionCoefficient = 0;

    _nMeasuredStars = 0;
    // The various _npar... are initialized in assignIndices.
    // Although there is no reason to adress them before one might be tempted by
    // evaluating a Chi2 rightaway, .. which uses these counts, so:
    assignIndices("");
}

#define NPAR_PM 2

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstrometryFit::transformFittedStar(const FittedStar &fittedStar, const Gtransfo *sky2TP,
                                         const Point &refractionVector, double refractionCoeff,
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
static void tweakAstromMeasurementErrors(FatPoint &P, const MeasuredStar &Ms, double error) {
    static bool called = false;
    static double increment = 0;
    if (!called) {
        increment = sqr(error);  // was in Preferences
        called = true;
    }
    P.vx += increment;
    P.vy += increment;
}

// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstrometryFit::LSDerivatives1(const CcdImage &ccdImage, TripletList &tList, Eigen::VectorXd &rhs,
                                   const MeasuredStarList *msList) const {
    /***************************************************************************/
    /**  Changes in this routine should be reflected into accumulateStatImage  */
    /***************************************************************************/
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
    vector<unsigned> indices(npar_tot, -1);
    if (_fittingDistortions) mapping->setMappingIndices(indices);

    // proper motion stuff
    double mjd = ccdImage.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = ccdImage.getRefractionVector();
    // transformation from sky to TP
    const Gtransfo *sky2TP = _astrometryModel->getSky2TP(ccdImage);
    // reserve matrices once for all measurements
    GtransfoLin dypdy;
    // the shape of h (et al) is required this way in order to be able to
    // separate derivatives along x and y as vectors.
    Eigen::MatrixX2d h(npar_tot, 2), halpha(npar_tot, 2), hw(npar_tot, 2);
    Eigen::Matrix2d transW(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::VectorXd grad(npar_tot);
    // current position in the Jacobian
    unsigned kTriplets = tList.getNextFreeIndex();
    const MeasuredStarList &catalog = (msList) ? *msList : ccdImage.getCatalogForFit();

    for (auto &i : catalog) {
        const MeasuredStar &ms = *i;
        if (!ms.isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = ms;
        tweakAstromMeasurementErrors(inPos, ms, _posError);
        h.setZero();  // we cannot be sure that all entries will be overwritten.
        FatPoint outPos;
        // should *not* fill h if whatToFit excludes mapping parameters.
        if (_fittingDistortions)
            mapping->computeTransformAndDerivatives(inPos, outPos, h);
        else
            mapping->transformPosAndErrors(inPos, outPos);

        unsigned ipar = npar_mapping;
        double det = outPos.vx * outPos.vy - sqr(outPos.vxy);
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
            h(npar_mapping, 0) = -dypdy.A11();
            h(npar_mapping + 1, 0) = -dypdy.A12();
            h(npar_mapping, 1) = -dypdy.A21();
            h(npar_mapping + 1, 1) = -dypdy.A22();
            indices[npar_mapping] = fs->getIndexInMatrix();
            indices.at(npar_mapping + 1) = fs->getIndexInMatrix() + 1;
            ipar += npar_pos;
        }
        /* only consider proper motions of objects allowed to move,
        unless the fit is going to be degenerate */
        if (_fittingPM && fs->mightMove) {
            h(ipar, 0) = -mjd;  // Sign unchecked but consistent with above
            h(ipar + 1, 1) = -mjd;
            indices[ipar] = fs->getIndexInMatrix() + 2;
            indices[ipar + 1] = fs->getIndexInMatrix() + 3;
            ipar += npar_pm;
        }
        if (_fittingRefrac) {
            /* if the definition of color changes, it has to remain
               consistent with transformFittedStar */
            double color = fs->color - _referenceColor;
            // sign checked
            h(ipar, 0) = -refractionVector.x * color;
            h(ipar, 1) = -refractionVector.y * color;
            indices[ipar] = _refracPosInMatrix;
            ipar += 1;
        }

        // We can now compute the residual
        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);

        // do not write grad = h*transW*res to avoid
        // dynamic allocation of a temporary
        halpha = h * alpha;
        hw = h * transW;
        grad = hw * res;
        // now feed in triplets and rhs
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar) {
            for (unsigned ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            rhs(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }                    // end loop on measurements
    tList.setNextFreeIndex(kTriplets);
}

void AstrometryFit::LSDerivatives2(const FittedStarList &fsl, TripletList &tList,
                                   Eigen::VectorXd &rhs) const {
    /* We compute here the derivatives of the terms involving fitted
       stars and reference stars. They only provide contributions if we
       are fitting positions: */
    if (!_fittingPos) return;
    /* the other case where the accumulation of derivatives stops
       here is when there are no RefStars */
    if (_associations.refStarList.size() == 0) return;
    Eigen::Matrix2d w(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::Matrix2d h(2, 2), halpha(2, 2), hw(2, 2);
    GtransfoLin der;
    Eigen::Vector2d res, grad;
    unsigned indices[2 + NPAR_PM];
    unsigned kTriplets = tList.getNextFreeIndex();
    /* We cannot use the spherical coordinates directly to evaluate
       Euclidean distances, we have to use a projector on some plane in
       order to express least squares. Not projecting could lead to a
       disaster around the poles or across alpha=0.  So we need a
       projector. We construct a projector and will change its
       projection point at every object */
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i : fsl) {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        proj.computeDerivative(fs, der, 1e-4);
        // sign checked. TODO check that the off-diagonal terms are OK.
        h(0, 0) = -der.A11();
        h(1, 0) = -der.A12();
        h(0, 1) = -der.A21();
        h(1, 1) = -der.A22();
        // TO DO : account for proper motions.
        double det = rsProj.vx * rsProj.vy - sqr(rsProj.vxy);
        if (rsProj.vx <= 0 || rsProj.vy <= 0 || det <= 0) {
            LOGLS_WARN(_log, "RefStar error matrix not positive definite for:  " << *rs);
            continue;
        }
        w(0, 0) = rsProj.vy / det;
        w(0, 1) = w(1, 0) = -rsProj.vxy / det;
        w(1, 1) = rsProj.vx / det;
        // compute alpha, a triangular square root
        // of w (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(w(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = w(0, 1) / alpha(0, 0);
        alpha(1, 1) = 1. / sqrt(det * w(0, 0));
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
        halpha = h * alpha;
        // grad = h*w*res
        hw = h * w;
        grad = hw * res;
        // now feed in triplets and rhs
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar) {
            for (unsigned ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            rhs(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }
    tList.setNextFreeIndex(kTriplets);
}

//! this routine computes the derivatives of all LS terms, including the ones that refer to references stars,
//! if any
void AstrometryFit::LSDerivatives(TripletList &tList, Eigen::VectorXd &rhs) const {
    auto L = _associations.getCcdImageList();
    for (auto const &im : L) {
        LSDerivatives1(*im, tList, rhs);
    }
    LSDerivatives2(_associations.fittedStarList, tList, rhs);
}

// This is almost a selection of lines of LSDerivatives1(CcdImage ...)
/* This routine (and the following one) is template because it is used
both with its first argument as "const CCdImage &" and "CcdImage &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ImType, class Accum>
void AstrometryFit::accumulateStatImage(ImType &image, Accum &accu) const {
    /**********************************************************************/
    /**  Changes in this routine should be reflected into LSDerivatives1  */
    /**********************************************************************/
    /* Setup */
    // 1 : get the Mapping's
    const Mapping *mapping = _astrometryModel->getMapping(image);
    // proper motion stuff
    double mjd = image.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = image.getRefractionVector();
    // transformation from sky to TP
    const Gtransfo *sky2TP = _astrometryModel->getSky2TP(image);
    // reserve matrix once for all measurements
    Eigen::Matrix2Xd transW(2, 2);

    auto &catalog = image.getCatalogForFit();
    for (auto const &ms : catalog) {
        if (!ms->isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = *ms;
        tweakAstromMeasurementErrors(inPos, *ms, _posError);

        FatPoint outPos;
        // should *not* fill h if whatToFit excludes mapping parameters.
        mapping->transformPosAndErrors(inPos, outPos);
        double det = outPos.vx * outPos.vy - sqr(outPos.vxy);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, " Inconsistent measurement errors :drop measurement at "
                                     << Point(*ms) << " in image " << image.getName());
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

        accu.addEntry(chi2Val, 2, ms);
    }  // end of loop on measurements
}

//! for a list of images.
template <class ListType, class Accum>
void AstrometryFit::accumulateStatImageList(ListType &list, Accum &accum) const {
    for (auto &im : list) {
        accumulateStatImage(*im, accum);
    }
}

template <class Accum>
void AstrometryFit::accumulateStatRefStars(Accum &accum) const {
    /* If you wonder why we project here, read comments in
       AstrometryFit::LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) */
    FittedStarList &fsl = _associations.fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &fs : fsl) {
        const RefStar *rs = fs->getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(*fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        // TO DO : account for proper motions.
        double rx = rsProj.x;  // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx * rsProj.vy - sqr(rsProj.vxy);
        double wxx = rsProj.vy / det;
        double wyy = rsProj.vx / det;
        double wxy = -rsProj.vxy / det;
        accum.addEntry(wxx * sqr(rx) + 2 * wxy * rx * ry + wyy * sqr(ry), 2, fs);
    }
}

Chi2Statistic AstrometryFit::computeChi2() const {
    Chi2Statistic chi2;
    accumulateStatImageList(_associations.getCcdImageList(), chi2);
    // now ref stars:
    accumulateStatRefStars(chi2);
    // so far, ndof contains the number of squares.
    // So, subtract here the number of parameters.
    chi2.ndof -= _nParTot;
    return chi2;
}

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void AstrometryFit::setMeasuredStarIndices(const MeasuredStar &ms, std::vector<unsigned> &indices) const {
    if (_fittingDistortions) {
        const Mapping *mapping = _astrometryModel->getMapping(ms.getCcdImage());
        mapping->setMappingIndices(indices);
    }
    auto fs = ms.getFittedStar();
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

//! contributions to derivatives of (presumambly) outlier terms. No discarding done.
void AstrometryFit::outliersContributions(MeasuredStarList &msOutliers, FittedStarList &fOutliers,
                                          TripletList &tList, Eigen::VectorXd &grad) {
    // contributions from measurement terms:
    for (auto const &i : msOutliers) {
        MeasuredStarList tmp;
        tmp.push_back(i);
        const CcdImage &ccd = i->getCcdImage();
        LSDerivatives1(ccd, tList, grad, &tmp);
    }
    LSDerivatives2(fOutliers, tList, grad);
}

unsigned AstrometryFit::removeOutliers(double nSigmaCut, const std::string &measOrRef) {
    MeasuredStarList msOutliers;
    FittedStarList fsOutliers;
    unsigned n = findOutliers(nSigmaCut, msOutliers, fsOutliers, measOrRef);
    removeMeasOutliers(msOutliers);
    removeRefOutliers(fsOutliers);
    return n;
}

unsigned AstrometryFit::findOutliers(double nSigmaCut, MeasuredStarList &msOutliers,
                                     FittedStarList &fsOutliers, const std::string &measOrRef) const {
    bool searchMeas = (measOrRef.find("Meas") != std::string::npos);
    bool searchRef = (measOrRef.find("Ref") != std::string::npos);

    // collect chi2 contributions
    Chi2List chi2List;
    chi2List.reserve(_nMeasuredStars + _associations.refStarList.size());
    // contributions from measurement terms:
    if (searchMeas) accumulateStatImageList(_associations.ccdImageList, chi2List);
    // and from reference terms
    if (searchRef) accumulateStatRefStars(chi2List);

    // compute some statistics
    size_t nval = chi2List.size();
    if (nval == 0) return 0;
    sort(chi2List.begin(), chi2List.end());
    double median = (nval & 1) ? chi2List[nval / 2].chi2
                               : 0.5 * (chi2List[nval / 2 - 1].chi2 + chi2List[nval / 2].chi2);
    auto averageAndSigma = chi2List.computeAverageAndSigma();
    LOGLS_DEBUG(_log, "RemoveOutliers chi2 stat: mean/median/sigma " << averageAndSigma.first << '/' << median
                                                                     << '/' << averageAndSigma.second);
    double cut = averageAndSigma.first + nSigmaCut * averageAndSigma.second;
    /* For each of the parameters, we will not remove more than 1
       measurement that contributes to constraining it. Keep track using
       of what we are touching using an integer vector. This is the
       trick that Marc Betoule came up to for outlier removals in "star
       flats" fits. */
    Eigen::VectorXi affectedParams(_nParTot);
    affectedParams.setZero();

    unsigned nOutliers = 0;  // returned to the caller
    // start from the strongest outliers.
    for (auto chi2 = chi2List.rbegin(); chi2 != chi2List.rend(); ++chi2) {
        if (chi2->chi2 < cut) break;  // because the array is sorted.
        vector<unsigned> indices;
        indices.reserve(100);  // just there to limit reallocations.
        /* now, we want to get the indices of the parameters this chi2
        term depends on. We have to figure out which kind of term it
         is; we use for that the type of the star attached to the Chi2Star. */
        auto ms = std::dynamic_pointer_cast<MeasuredStar>(chi2->star);
        std::shared_ptr<FittedStar> fs;
        if (!ms)  // it is reference term.
        {
            fs = std::dynamic_pointer_cast<FittedStar>(chi2->star);
            indices.push_back(fs->getIndexInMatrix());
            indices.push_back(fs->getIndexInMatrix() + 1);  // probably useless
            /* One might think it would be useful to account for PM
               parameters here, but it is just useless */
        } else  // it is a measurement term.
        {
            setMeasuredStarIndices(*ms, indices);
        }

        /* Find out if we already discarded a stronger outlier
        constraining some parameter this one constrains as well. If
         yes, we keep this one, because this stronger outlier could be
         causing the large chi2 we have in hand.  */
        bool drop_it = true;
        for (auto const &i : indices)
            if (affectedParams(i) != 0) drop_it = false;

        if (drop_it)  // store the outlier in one of the lists:
        {
            if (ms)  // measurement term
                msOutliers.push_back(ms);
            else  // ref term
                fsOutliers.push_back(fs);
            // mark the parameters as directly changed when we discard this chi2 term.
            for (auto const &i : indices) affectedParams(i)++;
            nOutliers++;
        }
    }  // end loop on measurements/references
    LOGLS_INFO(_log, "findOutliers: found " << msOutliers.size() << " meas outliers and " << fsOutliers.size()
                                            << " ref outliers ");

    return nOutliers;
}

void AstrometryFit::removeMeasOutliers(MeasuredStarList &outliers) {
    for (auto &i : outliers) {
        MeasuredStar &ms = *i;
        auto fs = std::const_pointer_cast<FittedStar>(ms.getFittedStar());
        ms.setValid(false);
        fs->getMeasurementCount()--;  // could be put in setValid
    }
}

void AstrometryFit::removeRefOutliers(FittedStarList &outliers) {
    for (auto &i : outliers) {
        FittedStar &fs = *i;
        fs.setRefStar(nullptr);
    }
}

void AstrometryFit::assignIndices(const std::string &whatToFit) {
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: Now fitting " << whatToFit);
    _fittingDistortions = (_whatToFit.find("Distortions") != string::npos);
    _fittingPos = (_whatToFit.find("Positions") != string::npos);
    _fittingRefrac = (_whatToFit.find("Refrac") != string::npos);
    if (_sigCol == 0 && _fittingRefrac) {
        LOGLS_WARN(_log,
                   "Cannot fit refraction coefficients without a color lever arm. Ignoring refraction.");
        _fittingRefrac = false;
    }
    _fittingPM = (_whatToFit.find("PM") != string::npos);
    // When entering here, we assume that whatToFit has already been interpreted.

    _nParDistortions = 0;
    if (_fittingDistortions) _nParDistortions = _astrometryModel->assignIndices(0, _whatToFit);
    unsigned ipar = _nParDistortions;

    if (_fittingPos) {
        FittedStarList &fsl = _associations.fittedStarList;
        for (auto const &i : fsl) {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in GetMeasuredStarIndices
            fs.setIndexInMatrix(ipar);
            ipar += 2;
            if ((_fittingPM)&fs.mightMove) ipar += NPAR_PM;
        }
    }
    _nParPositions = ipar - _nParDistortions;
    if (_fittingRefrac) {
        _refracPosInMatrix = ipar;
        ipar += _nRefrac;
    }
    _nParTot = ipar;
}

void AstrometryFit::offsetParams(const Eigen::VectorXd &delta) {
    if (delta.size() != _nParTot)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "AstrometryFit::offsetParams : the provided vector length is not compatible with "
                          "the current whatToFit setting");
    if (_fittingDistortions) _astrometryModel->offsetParams(delta);

    if (_fittingPos) {
        FittedStarList &fsl = _associations.fittedStarList;
        for (auto const &i : fsl) {
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
static void write_sparse_matrix_in_fits(const SpMat &mat, const string &fitsName) {
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

static void write_vect_in_fits(const Eigen::VectorXd &vectorXd, const string &fitsName) {
    Vect v(vectorXd.size());
    for (int k = 0; k < vectorXd.size(); ++k) v(k) = V(k);
    Mat(v).writeFits(fitsName);
}

#endif

int AstrometryFit::minimize(const std::string &whatToFit, const double nSigRejCut) {
    assignIndices(whatToFit);

    // return code can take 3 values :
    // 0 : fit has converged - no more outliers
    // 1 : still some ouliers but chi2 increases
    // 2 : factorization failed
    int returnCode = 0;

    // TODO : write a guesser for the number of triplets
    unsigned nTrip = (_LastNTrip) ? _LastNTrip : 1e6;
    TripletList tList(nTrip);
    Eigen::VectorXd grad(_nParTot);
    grad.setZero();

    // Fill the triplets
    LSDerivatives(tList, grad);
    _LastNTrip = tList.size();

    LOGLS_DEBUG(_log, "End of triplet filling, ntrip = " << tList.size());

    SpMat hessian;
    {
        SpMat jacobian(_nParTot, tList.getNextFreeIndex());
        jacobian.setFromTriplets(tList.begin(), tList.end());
        // release memory shrink_to_fit is C++11
        tList.clear();  // tList.shrink_to_fit();
        hessian = jacobian * jacobian.transpose();
    }  // release the Jacobian

    LOGLS_DEBUG(_log, "Starting factorization, hessian: dim="
                              << hessian.rows() << " nnz=" << hessian.nonZeros()
                              << " filling-frac = " << hessian.nonZeros() / sqr(hessian.rows()));

    CholmodSimplicialLDLT2<SpMat> chol(hessian);
    if (chol.info() != Eigen::Success) {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return 2;
    }

    unsigned tot_outliers = 0;
    double old_chi2 = computeChi2().chi2;

    while (true) {
        Eigen::VectorXd delta = chol.solve(grad);
        offsetParams(delta);
        Chi2Statistic current_chi2(computeChi2());
        LOGLS_DEBUG(_log, current_chi2);
        if (current_chi2.chi2 > old_chi2) {
            LOGL_WARN(_log, "chi2 went up, exiting outlier rejection loop");
            returnCode = 1;
            break;
        }
        old_chi2 = current_chi2.chi2;

        if (nSigRejCut == 0) break;
        MeasuredStarList moutliers;
        FittedStarList foutliers;
        int n_outliers = findOutliers(nSigRejCut, moutliers, foutliers);
        tot_outliers += n_outliers;
        if (n_outliers == 0) break;
        TripletList tList(1000);  // initial allocation size.
        grad.setZero();           // recycle the gradient
        // compute the contributions of outliers to derivatives
        outliersContributions(moutliers, foutliers, tList, grad);
        // actually discard them
        removeMeasOutliers(moutliers);
        removeRefOutliers(foutliers);
        // convert triplet list to eigen internal format
        SpMat h(_nParTot, tList.getNextFreeIndex());
        h.setFromTriplets(tList.begin(), tList.end());
        int update_status = chol.update(h, false /* means downdate */);
        LOGLS_DEBUG(_log, "cholmod update_status " << update_status);
        /* The contribution of outliers to the gradient is the opposite
        of the contribution of all other terms, because they add up
         to 0 */
        grad *= -1;
    }

    LOGLS_INFO(_log, "Total number of outliers " << tot_outliers);

    return returnCode;
}

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
        TripletList tList(10000);
        Eigen::VectorXd rhs(_nParTot);
        rhs.setZero();
        LSDerivatives(tList, rhs);
        SpMat jacobian(_nParTot, tList.getNextFreeIndex());
        jacobian.setFromTriplets(tList.begin(), tList.end());
        SpMat hessian = jacobian * jacobian.transpose();
#ifdef STORAGE
        char name[24];
        sprintf(name, "h%d.fits", k);
        write_sparse_matrix_in_fits(hessian, name);
        sprintf(name, "g%d.fits", k);
        write_vect_in_fits(rhs, name);
#endif
        LOGLS_DEBUG(_log, "npar : " << _nParTot << ' ' << _nParDistortions);
    }
}

void AstrometryFit::makeResTuple(const std::string &tupleName) const {
    /* cook-up 2 different file names by inserting something just before
       the dot (if any), and within the actual file name. */
    size_t dot = tupleName.rfind('.');
    size_t slash = tupleName.rfind('/');
    if (dot == string::npos || (slash != string::npos && dot < slash)) dot = tupleName.size();
    std::string meas_tuple(tupleName);
    meas_tuple.insert(dot, "-meas");
    makeMeasResTuple(meas_tuple);
    std::string ref_tuple(tupleName);
    ref_tuple.insert(dot, "-ref");
    makeRefResTuple(ref_tuple);
}

void AstrometryFit::makeMeasResTuple(const std::string &tupleName) const {
    std::ofstream tuple(tupleName.c_str());
    tuple << "#xccd: coordinate in CCD" << endl
          << "#yccd: " << endl
          << "#rx:   residual in degrees in TP" << endl
          << "#ry:" << endl
          << "#xtp: transformed coordinate in TP " << endl
          << "#ytp:" << endl
          << "#mag: rough mag" << endl
          << "#jd: Julian date of the measurement" << endl
          << "#rvx: transformed measurement uncertainty " << endl
          << "#rvy:" << endl
          << "#rvxy:" << endl
          << "#color : " << endl
          << "#fsindex: some unique index of the object" << endl
          << "#ra: pos of fitted star" << endl
          << "#dec: pos of fitted star" << endl
          << "#chi2: contribution to Chi2 (2D dofs)" << endl
          << "#nm: number of measurements of this FittedStar" << endl
          << "#chip: chip number" << endl
          << "#visit: visit id" << endl
          << "#end" << endl;
    const CcdImageList &L = _associations.getCcdImageList();
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
            double det = tpPos.vx * tpPos.vy - sqr(tpPos.vxy);
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
                  << im.getVisit() << endl;
        }  // loop on measurements in image
    }      // loop on images
}

void AstrometryFit::makeRefResTuple(const std::string &tupleName) const {
    std::ofstream tuple(tupleName.c_str());
    tuple << "#ra: coordinates of FittedStar" << endl
          << "#dec: " << endl
          << "#rx:   residual in degrees in TP" << endl
          << "#ry:" << endl
          << "#mag: mag" << endl
          << "#rvx: transformed measurement uncertainty " << endl
          << "#rvy:" << endl
          << "#rvxy:" << endl
          << "#color : " << endl
          << "#fsindex: some unique index of the object" << endl
          << "#chi2: contribution to Chi2 (2D dofs)" << endl
          << "#nm: number of measurements of this FittedStar" << endl
          << "#end" << endl;
    // The following loop is heavily inspired from AstrometryFit::computeChi2()
    const FittedStarList &fsl = _associations.fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i : fsl) {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        double rx = rsProj.x;  // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx * rsProj.vy - sqr(rsProj.vxy);
        double wxx = rsProj.vy / det;
        double wyy = rsProj.vx / det;
        double wxy = -rsProj.vxy / det;
        double chi2 = wxx * sqr(rx) + 2 * wxy * rx * ry + wyy * sqr(ry);
        tuple << std::setprecision(9);
        tuple << fs.x << ' ' << fs.y << ' ' << rx << ' ' << ry << ' ' << fs.getMag() << ' ' << rsProj.vx
              << ' ' << rsProj.vy << ' ' << rsProj.vxy << ' ' << fs.color << ' ' << fs.getIndexInMatrix()
              << ' ' << chi2 << ' ' << fs.getMeasurementCount() << endl;
    }  // loop on FittedStars
}
}  // namespace jointcal
}  // namespace lsst
