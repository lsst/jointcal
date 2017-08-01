#include <vector>

#include "lsst/log/Log.h"

#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.Fitter");
}

namespace lsst {
namespace jointcal {

Chi2Statistic FitterBase::computeChi2() const {
    Chi2Statistic chi2;
    accumulateStatImageList(_associations->getCcdImageList(), chi2);
    accumulateStatRefStars(chi2);
    // chi2.ndof contains the number of squares.
    // So subtract the number of parameters.
    chi2.ndof -= _nParTot;
    return chi2;
}

unsigned FitterBase::findOutliers(double nSigmaCut, MeasuredStarList &msOutliers,
                                  FittedStarList &fsOutliers) const {
    // collect chi2 contributions
    Chi2List chi2List;
    chi2List.reserve(_nMeasuredStars + _associations->refStarList.size());
    // contributions from measurement terms:
    accumulateStatImageList(_associations->ccdImageList, chi2List);
    // and from reference terms
    accumulateStatRefStars(chi2List);

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
        std::vector<unsigned> indices;
        /* now, we want to get the indices of the parameters this chi2
        term depends on. We have to figure out which kind of term it
         is; we use for that the type of the star attached to the Chi2Star. */
        auto ms = std::dynamic_pointer_cast<MeasuredStar>(chi2->star);
        std::shared_ptr<FittedStar> fs;
        if (!ms) {
            // it is reference term.
            fs = std::dynamic_pointer_cast<FittedStar>(chi2->star);
            indices.push_back(fs->getIndexInMatrix());
            indices.push_back(fs->getIndexInMatrix() + 1);  // probably useless
            /* One might think it would be useful to account for PM
               parameters here, but it is just useless */
        } else {  // it is a measurement term.
            getIndicesOfMeasuredStar(*ms, indices);
        }

        /* Find out if we already discarded a stronger outlier
        constraining some parameter this one constrains as well. If
         yes, we keep this one, because this stronger outlier could be
         causing the large chi2 we have in hand.  */
        bool drop_it = true;
        for (auto const &i : indices) {
            if (affectedParams(i) != 0) {
                drop_it = false;
            }
        }

        if (drop_it)  // store the outlier in one of the lists:
        {
            if (ms) {
                // measurement term
                msOutliers.push_back(ms);
            } else {
                // ref term
                fsOutliers.push_back(fs);
            }
            // mark the parameters as directly changed when we discard this chi2 term.
            for (auto const &i : indices) {
                affectedParams(i)++;
            }
            nOutliers++;
        }
    }  // end loop on measurements/references
    LOGLS_INFO(_log, "findOutliers: found " << msOutliers.size() << " meas outliers and " << fsOutliers.size()
                                            << " ref outliers ");

    return nOutliers;
}

MinimizeResult FitterBase::minimize(std::string const &whatToFit, double nSigmaCut) {
    assignIndices(whatToFit);

    MinimizeResult returnCode = MinimizeResult::Converged;

    // TODO : write a guesser for the number of triplets
    unsigned nTrip = (_lastNTrip) ? _lastNTrip : 1e6;
    TripletList tripletList(nTrip);
    Eigen::VectorXd grad(_nParTot);
    grad.setZero();

    // Fill the triplets
    leastSquareDerivatives(tripletList, grad);
    _lastNTrip = tripletList.size();

    LOGLS_DEBUG(_log, "End of triplet filling, ntrip = " << tripletList.size());

    SpMat hessian;
    {
        SpMat jacobian(_nParTot, tripletList.getNextFreeIndex());
        jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
        // release memory shrink_to_fit is C++11
        tripletList.clear();  // tripletList.shrink_to_fit();
        hessian = jacobian * jacobian.transpose();
    }  // release the Jacobian

    LOGLS_DEBUG(_log, "Starting factorization, hessian: dim="
                              << hessian.rows() << " non-zeros=" << hessian.nonZeros()
                              << " filling-frac = " << hessian.nonZeros() / std::pow(hessian.rows(), 2));

    CholmodSimplicialLDLT2<SpMat> chol(hessian);
    if (chol.info() != Eigen::Success) {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return MinimizeResult::Failed;
    }

    unsigned totalOutliers = 0;
    double oldChi2 = computeChi2().chi2;

    while (true) {
        Eigen::VectorXd delta = chol.solve(grad);
        offsetParams(delta);
        Chi2Statistic currentChi2(computeChi2());
        LOGLS_DEBUG(_log, currentChi2);
        if (currentChi2.chi2 > oldChi2) {
            LOGL_WARN(_log, "chi2 went up, skipping outlier rejection loop");
            returnCode = MinimizeResult::Chi2Increased;
            break;
        }
        oldChi2 = currentChi2.chi2;

        if (nSigmaCut == 0) break;  // no rejection step to perform
        MeasuredStarList msOutliers;
        FittedStarList fsOutliers;
        int nOutliers = findOutliers(nSigmaCut, msOutliers, fsOutliers);
        totalOutliers += nOutliers;
        if (nOutliers == 0) break;
        TripletList tripletList(nOutliers);
        grad.setZero();  // recycle the gradient
        // compute the contributions of outliers to derivatives
        outliersContributions(msOutliers, fsOutliers, tripletList, grad);
        // Remove significant outliers
        removeMeasOutliers(msOutliers);
        removeRefOutliers(fsOutliers);
        // convert triplet list to eigen internal format
        SpMat H(_nParTot, tripletList.getNextFreeIndex());
        H.setFromTriplets(tripletList.begin(), tripletList.end());
        int update_status = chol.update(H, false /* means downdate */);
        LOGLS_DEBUG(_log, "cholmod update_status " << update_status);
        // The contribution of outliers to the gradient is the opposite
        // of the contribution of all other terms, because they add up to 0
        grad *= -1;
    }

    LOGLS_INFO(_log, "Total number of outliers " << totalOutliers);
    return returnCode;
}

void FitterBase::outliersContributions(MeasuredStarList &msOutliers, FittedStarList &fsOutliers,
                                       TripletList &tripletList, Eigen::VectorXd &grad) {
    for (auto &outlier : msOutliers) {
        MeasuredStarList tmp;
        tmp.push_back(outlier);
        const CcdImage &ccdImage = outlier->getCcdImage();
        leastSquareDerivativesMeasurement(ccdImage, tripletList, grad, &tmp);
    }
    leastSquareDerivativesReference(fsOutliers, tripletList, grad);
}

void FitterBase::removeMeasOutliers(MeasuredStarList &outliers) {
    for (auto &measuredStar : outliers) {
        auto fittedStar = std::const_pointer_cast<FittedStar>(measuredStar->getFittedStar());
        measuredStar->setValid(false);
        fittedStar->getMeasurementCount()--;  // could be put in setValid
    }
}

void FitterBase::removeRefOutliers(FittedStarList &outliers) {
    for (auto &fittedStar : outliers) {
        fittedStar->setRefStar(nullptr);
    }
}

void FitterBase::leastSquareDerivatives(TripletList &tripletList, Eigen::VectorXd &grad) const {
    auto ccdImageList = _associations->getCcdImageList();
    for (auto const &ccdImage : ccdImageList) {
        leastSquareDerivativesMeasurement(*ccdImage, tripletList, grad);
    }
    leastSquareDerivativesReference(_associations->fittedStarList, tripletList, grad);
}

}  // namespace jointcal
}  // namespace lsst
