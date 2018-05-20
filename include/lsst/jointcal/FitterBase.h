// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FITTER_BASE_H
#define LSST_JOINTCAL_FITTER_BASE_H

#include "lsst/log/Log.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Tripletlist.h"

namespace lsst {
namespace jointcal {

/// Return value of minimize()
enum class MinimizeResult {
    Converged,      // fit has converged - no more outliers
    Chi2Increased,  // still some ouliers but chi2 increases
    Failed          // factorization failed
};

/**
 * Base class for fitters.
 *
 * Implements minimize and findOutliers. Chi2, residual, derivative, etc. calculations must be
 * implemented in the child class via the virtual methods.
 */
class FitterBase {
public:
    explicit FitterBase(std::shared_ptr<Associations> associations)
            : _associations(associations), _whatToFit(""), _lastNTrip(0), _nParTot(0), _nMeasuredStars(0) {}

    /// No copy or move: there is only ever one fitter of a given type.
    FitterBase(FitterBase const &) = delete;
    FitterBase(FitterBase &&) = delete;
    FitterBase &operator=(FitterBase const &) = delete;
    FitterBase &operator=(FitterBase &&) = delete;

    /**
     * Does a 1 step minimization, assuming a linear model.
     *
     * This is a complete Newton Raphson step. Compute first and second
     * derivatives, solve for the step and apply it, without a line search.
     *
     * It calls assignIndices, leastSquareDerivatives, solves the linear system and calls
     * offsetParams, then removes outliers in a loop if requested.
     * Relies on sparse linear algebra.
     *
     * @param[in]  whatToFit  See child method assignIndices for valid string values.
     * @param[in]  nSigmaCut  How many sigma to reject outliers at. Outlier
     *                        rejection ignored for nSigmaCut=0.
     * @param[in]  doRankUpdate  Use CholmodSimplicialLDLT2.update() to do a fast rank update after outlier
     *                           removal; otherwise do a slower full recomputation of the matrix.
     *                           Only matters if nSigmaCut != 0.
     *
     * @return  Return code describing success/failure of fit.
     *
     * @note   When fitting one parameter set by itself (e.g. "Model"), the system is purely linear,
     *         which should result in the optimal chi2 after a single step. This can
     *         be used to debug the fitter by fitting that parameter set twice in a row:
     *         the second run with the same "whatToFit" will produce no change in
     *         the fitted parameters, if the calculations and indices are defined correctly.
     */
    MinimizeResult minimize(std::string const &whatToFit, double const nSigmaCut = 0,
                            bool const doRankUpdate = true);

    /**
     * Returns the chi2 for the current state.
     */
    Chi2Statistic computeChi2() const;

    /**
     * Evaluates the chI^2 derivatives (Jacobian and gradient) for the current whatToFit setting.
     *
     * The Jacobian is given as triplets in a sparse matrix, the gradient as a dense vector.
     * The parameters which vary, and their indices, are to be set using  assignIndices.
     *
     * @param      tripletList  tripletList of (row,col,value) representing the Jacobian of the chi2.
     * @param      grad         The gradient of the chi2.
     */
    void leastSquareDerivatives(TripletList &tripletList, Eigen::VectorXd &grad) const;

    /**
     * Offset the parameters by the requested quantities. The used parameter
     * layout is the one from the last call to assignIndices or minimize(). There
     * is no easy way to check that the current setting of whatToFit and the
     * provided Delta vector are compatible: we can only test the size.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    /**
     * Set parameters to fit and assign indices in the big matrix.
     *
     * @param[in]  whatToFit  See child class documentation for valid string values.
     */
    virtual void assignIndices(std::string const &whatToFit) = 0;

    /**
     * Save the full chi2 term per star that was used in the minimization, for debugging.
     *
     * Saves results to text files "baseName-meas.csv" and "baseName-ref.csv" for the
     * MeasuredStar and RefStar contributions, respectively.
     * This method is mostly useful for debugging: we will probably want to create a better persistence
     * system for jointcal's internal representations in the future (see DM-12446).
     */
    virtual void saveChi2Contributions(std::string const &baseName) const;

    /// Save a CSV file containing residuals of measurement terms.
    virtual void saveChi2MeasContributions(std::string const &baseName) const = 0;

    /// Save a CSV file containing residuals of reference terms.
    virtual void saveChi2RefContributions(std::string const &baseName) const = 0;

protected:
    std::shared_ptr<Associations> _associations;
    std::string _whatToFit;

    int _lastNTrip;  // last triplet count, used to speed up allocation
    unsigned int _nParTot;
    unsigned _nMeasuredStars;

    // lsst.logging instance, to be created by subclass so that messages have consistent name while fitting.
    LOG_LOGGER _log;

    /**
     * Find Measurements and references contributing more than a cut, computed as
     * @f[
     *     <chi2> + nSigmaCut + rms(chi2).
     * @f]
     * The outliers are NOT removed, and no refit is done.
     *
     * After returning from here, there are still measurements that contribute above the cut,
     * but their contribution should be evaluated after a refit before discarding them.
     *
     * @param[in]  nSigmaCut   Number of sigma to select on.
     * @param[out] msOutliers  list of MeasuredStar outliers to populate
     * @param[out] fsOutliers  list of FittedStar outliers to populate
     *
     * @return     Total number of outliers that were removed.
     */
    unsigned findOutliers(double nSigmaCut, MeasuredStarList &msOutliers, FittedStarList &fsOutliers) const;

    /**
     * Contributions to derivatives from (presumably) outlier terms. No
     * discarding done.
     */
    void outliersContributions(MeasuredStarList &msOutliers, FittedStarList &fsOutliers,
                               TripletList &tripletList, Eigen::VectorXd &grad);

    /// Remove measuredStar outliers from the fit. No Refit done.
    void removeMeasOutliers(MeasuredStarList &outliers);

    /// Remove refStar outliers from the fit. No Refit done.
    void removeRefOutliers(FittedStarList &outliers);

    /// Set the indices of a measured star from the full matrix, for outlier removal.
    virtual void getIndicesOfMeasuredStar(MeasuredStar const &measuredStar,
                                          std::vector<unsigned> &indices) const = 0;

    /// Compute the chi2 (per star or total, depending on which Chi2Accumulator is used) for measurements.
    virtual void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const = 0;

    /// Compute the chi2 (per star or total, depending on which Chi2Accumulator is used) for RefStars.
    virtual void accumulateStatRefStars(Chi2Accumulator &accum) const = 0;

    /**
     * Compute the derivatives of the measured stars and model for one CcdImage.
     *
     * The last argument will process a sub-list for outlier removal.
     */
    virtual void leastSquareDerivativesMeasurement(
            CcdImage const &ccdImage, TripletList &tripletList, Eigen::VectorXd &grad,
            MeasuredStarList const *measuredStarList = nullptr) const = 0;

    /// Compute the derivatives of the reference terms
    virtual void leastSquareDerivativesReference(FittedStarList const &fittedStarList,
                                                 TripletList &tripletList, Eigen::VectorXd &grad) const = 0;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_FITTER_BASE_H
