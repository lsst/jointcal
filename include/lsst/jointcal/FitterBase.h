// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FITTER_BASE_H
#define LSST_JOINTCAL_FITTER_BASE_H

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Tripletlist.h"

namespace lsst {
namespace jointcal {

/**
 * Base class for fitters.
 *
 * Implements minimize and findOutliers. Chi2, residual, derivative, etc. calculations must be
 * implemented in the child class via the virtual methods.
 */
class FitterBase {
public:
    FitterBase(Associations &associations) : _associations(associations), _lastNTrip(0), _nMeasuredStars(0) {
        // The various _npar are initialized in assignIndices, which must be called in the child constructor.
    }

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
     *
     * @return     Return code describing success of fit, can take 3 values:
     *             0 : fit has converged - no more outliers
     *             1 : still some ouliers but chi2 increased
     *             2 : factorization failed
     *
     * @note   When fitting one parameter set by itself (e.g. "Model"), the system is purely linear,
     *         which should result in the optimal chi2 after a single step. This can
     *         be used to debug the fitter by fitting that paramter set twice in a row:
     *         the second run with the same "whatToFit" will produce no change in
     *         the fitted parameters, if the calculations and indices are defined correctly.
     */
    int minimize(const std::string &whatToFit, double nSigmaCut = 0);

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
    virtual void offsetParams(const Eigen::VectorXd &delta) = 0;

    /**
     * Set parameters to fit and assign indices in the big matrix.
     *
     * @param[in]  whatToFit  See child class documentation for valid string values.
     */
    virtual void assignIndices(const std::string &whatToFit) = 0;

    /**
     * Save the full chi2 term per star that was used in the minimization.
     *
     * Produces both MeasuredStar and RefStar tuples.
     */
    virtual void saveResultTuples(const std::string &tupleName) const = 0;

protected:
    Associations &_associations;
    std::string _whatToFit;

    int _lastNTrip;  // last triplet count, used to speed up allocation
    unsigned int _nParTot;
    unsigned _nMeasuredStars;

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

    virtual void setMeasuredStarIndices(const MeasuredStar &measuredStar,
                                        std::vector<unsigned> &indices) const = 0;

    virtual void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const = 0;

    virtual void accumulateStatRefStars(Chi2Accumulator &accum) const = 0;

    /**
     * Compute the derivatives of the measured stars and model for one CcdImage.
     *
     * The last argument will process a sub-list for outlier removal.
     */
    virtual void leastSquareDerivativesMeasurement(
            const CcdImage &ccdImage, TripletList &tripletList, Eigen::VectorXd &grad,
            const MeasuredStarList *measuredStarList = nullptr) const = 0;

    /// Compute the derivatives of the reference terms
    virtual void leastSquareDerivativesReference(const FittedStarList &fittedStarList,
                                                 TripletList &tripletList, Eigen::VectorXd &grad) const = 0;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_FITTER_BASE_H
