// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <vector>
#include "Eigen/Core"

#include <boost/math/tools/minima.hpp>

#include "lsst/log/Log.h"

#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"

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
    LOGLS_DEBUG(_log, "findOutliers chi2 stat: mean/median/sigma " << averageAndSigma.first << '/' << median
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
        auto measuredStar = std::dynamic_pointer_cast<MeasuredStar>(chi2->star);
        std::shared_ptr<FittedStar> fittedStar;  // To add to fsOutliers if it is a reference outlier.
        if (measuredStar == nullptr) {
            // it is a reference outlier
            fittedStar = std::dynamic_pointer_cast<FittedStar>(chi2->star);
            if (fittedStar->getMeasurementCount() == 0) {
                LOGLS_WARN(_log, "FittedStar with no measuredStars found as an outlier: " << *fittedStar);
                continue;
            }
            // NOTE: Stars contribute twice to astrometry (x,y), but once to photometry (flux),
            // NOTE: but we only need to mark one index here because both will be removed with that star.
            indices.push_back(fittedStar->getIndexInMatrix());
            LOGLS_TRACE(_log, "Removing refStar " << *(fittedStar->getRefStar()) << " chi2: " << chi2->chi2);
            /* One might think it would be useful to account for PM
               parameters here, but it is just useless */
        } else {
            // it is a measurement outlier
            auto tempFittedStar = measuredStar->getFittedStar();
            if (tempFittedStar->getMeasurementCount() == 1 && tempFittedStar->getRefStar() == nullptr) {
                LOGLS_WARN(_log, "FittedStar with 1 measuredStar and no refStar found as an outlier: "
                                         << *tempFittedStar);
                continue;
            }
            getIndicesOfMeasuredStar(*measuredStar, indices);
            LOGLS_TRACE(_log, "Removing measStar " << *measuredStar << " chi2: " << chi2->chi2);
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
            if (measuredStar == nullptr) {
                // reference term
                fsOutliers.push_back(fittedStar);
            } else {
                // measurement term
                msOutliers.push_back(measuredStar);
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

namespace {
/// Return a Hessian matrix filled from tripletList of size nParTot x nParTot.
SparseMatrixD createHessian(int nParTot, TripletList const &tripletList) {
    SparseMatrixD jacobian(nParTot, tripletList.getNextFreeIndex());
    jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
    return jacobian * jacobian.transpose();
}

/// Write matrix and gradient to files built from dumpFile, and log their names.
void dumpMatrixAndGradient(SparseMatrixD const &matrix, Eigen::VectorXd const &grad,
                           std::string const &dumpFile, LOG_LOGGER _log) {
    std::string ext = ".txt";
    Eigen::MatrixXd matrixDense(matrix);
    std::string dumpMatrixPath = dumpFile + "-mat" + ext;
    std::ofstream matfile(dumpMatrixPath);
    matfile << matrixDense << std::endl;
    std::string dumpGradPath = dumpFile + "-grad" + ext;
    std::ofstream gradfile(dumpGradPath);
    gradfile << grad << std::endl;
    LOGLS_INFO(_log, "Dumped Hessian, gradient to: '" << dumpMatrixPath << "', '" << dumpGradPath << "'");
}
}  // namespace

MinimizeResult FitterBase::minimize(std::string const &whatToFit, double nSigmaCut, bool doRankUpdate,
                                    bool const doLineSearch, std::string const &dumpMatrixFile) {
    assignIndices(whatToFit);

    MinimizeResult returnCode = MinimizeResult::Converged;

    // TODO : write a guesser for the number of triplets
    unsigned nTrip = (_lastNTrip) ? _lastNTrip : 1e6;
    TripletList tripletList(nTrip);
    Eigen::VectorXd grad(_nParTot);
    grad.setZero();
    double scale = 1.0;

    // Fill the triplets
    leastSquareDerivatives(tripletList, grad);
    _lastNTrip = tripletList.size();

    LOGLS_DEBUG(_log, "End of triplet filling, ntrip = " << tripletList.size());

    SparseMatrixD hessian = createHessian(_nParTot, tripletList);
    tripletList.clear();  // we don't need it any more after we have the hessian.

    LOGLS_DEBUG(_log, "Starting factorization, hessian: dim="
                              << hessian.rows() << " non-zeros=" << hessian.nonZeros()
                              << " filling-frac = " << hessian.nonZeros() / std::pow(hessian.rows(), 2));

    if (dumpMatrixFile != "") {
        if (hessian.rows() * hessian.cols() > 2e8) {
            LOGLS_WARN(_log, "Hessian matrix is too big to dump to file, with rows, columns: "
                                     << hessian.rows() << ", " << hessian.cols());
        } else {
            dumpMatrixAndGradient(hessian, grad, dumpMatrixFile, _log);
        }
    }

    CholmodSimplicialLDLT2<SparseMatrixD> chol(hessian);
    if (chol.info() != Eigen::Success) {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return MinimizeResult::Failed;
    }

    unsigned totalMeasOutliers = 0;
    unsigned totalRefOutliers = 0;
    double oldChi2 = computeChi2().chi2;

    while (true) {
        Eigen::VectorXd delta = chol.solve(grad);
        if (doLineSearch) {
            scale = _lineSearch(delta);
        }
        offsetParams(scale * delta);
        Chi2Statistic currentChi2(computeChi2());
        LOGLS_DEBUG(_log, currentChi2);
        if (!isfinite(currentChi2.chi2)) {
            LOGL_ERROR(_log, "chi2 is not finite. Aborting outlier rejection.");
            returnCode = MinimizeResult::NonFinite;
            break;
        }
        if (currentChi2.chi2 > oldChi2 && totalMeasOutliers + totalRefOutliers != 0) {
            LOGL_WARN(_log, "chi2 went up, skipping outlier rejection loop");
            returnCode = MinimizeResult::Chi2Increased;
            break;
        }
        oldChi2 = currentChi2.chi2;

        if (nSigmaCut == 0) break;  // no rejection step to perform
        MeasuredStarList msOutliers;
        FittedStarList fsOutliers;
        // keep nOutliers so we don't have to sum msOutliers.size()+fsOutliers.size() twice below.
        int nOutliers = findOutliers(nSigmaCut, msOutliers, fsOutliers);
        totalMeasOutliers += msOutliers.size();
        totalRefOutliers += fsOutliers.size();
        if (nOutliers == 0) break;
        TripletList outlierTriplets(nOutliers);
        grad.setZero();  // recycle the gradient
        // compute the contributions of outliers to derivatives
        outliersContributions(msOutliers, fsOutliers, outlierTriplets, grad);
        // Remove significant outliers
        removeMeasOutliers(msOutliers);
        removeRefOutliers(fsOutliers);
        if (doRankUpdate) {
            // convert triplet list to eigen internal format
            SparseMatrixD H(_nParTot, outlierTriplets.getNextFreeIndex());
            H.setFromTriplets(outlierTriplets.begin(), outlierTriplets.end());
            chol.update(H, false /* means downdate */);
            // The contribution of outliers to the gradient is the opposite
            // of the contribution of all other terms, because they add up to 0
            grad *= -1;
        } else {
            // don't reuse tripletList because we want a new nextFreeIndex.
            TripletList nextTripletList(_lastNTrip);
            grad.setZero();
            // Rebuild the matrix and gradient
            leastSquareDerivatives(nextTripletList, grad);
            _lastNTrip = nextTripletList.size();
            LOGLS_DEBUG(_log, "Triplets recomputed, ntrip = " << nextTripletList.size());

            hessian = createHessian(_nParTot, nextTripletList);
            nextTripletList.clear();  // we don't need it any more after we have the hessian.

            LOGLS_DEBUG(_log,
                        "Restarting factorization, hessian: dim="
                                << hessian.rows() << " non-zeros=" << hessian.nonZeros()
                                << " filling-frac = " << hessian.nonZeros() / std::pow(hessian.rows(), 2));
            chol.compute(hessian);
            if (chol.info() != Eigen::Success) {
                LOGLS_ERROR(_log, "minimize: factorization failed ");
                return MinimizeResult::Failed;
            }
        }
    }

    // only print the outlier summary if outlier rejection was turned on.
    if (nSigmaCut != 0) {
        LOGLS_INFO(_log, "Number of outliers (Measured + Reference = Total): "
                                 << totalMeasOutliers << " + " << totalRefOutliers << " = "
                                 << totalMeasOutliers + totalRefOutliers);
    }
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
        auto fittedStar = measuredStar->getFittedStar();
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

void FitterBase::saveChi2Contributions(std::string const &baseName) const {
    std::string replaceStr = "{type}";
    auto pos = baseName.find(replaceStr);
    std::string measFilename(baseName);
    measFilename.replace(pos, replaceStr.size(), "-meas.csv");
    std::string refFilename(baseName);
    refFilename.replace(pos, replaceStr.size(), "-ref.csv");
    saveChi2MeasContributions(measFilename);
    saveChi2RefContributions(refFilename);
}

double FitterBase::_lineSearch(Eigen::VectorXd const &delta) {
    auto func = [this, &delta](double scale) {
        auto offset = scale * delta;
        offsetParams(offset);
        auto chi2 = computeChi2();
        // reset the system to where it was before offsetting.
        offsetParams(-offset);
        return chi2.chi2;
    };
    // The maximum theoretical precision is half the number of bits in the mantissa (see boost docs).
    auto bits = std::numeric_limits<double>::digits / 2;
    auto result = boost::math::tools::brent_find_minima(func, -1.0, 2.0, bits);
    LOGLS_DEBUG(_log, "Line search scale factor: " << result.first);
    return result.first;
}

}  // namespace jointcal
}  // namespace lsst
