#ifndef LSST_JOINTCAL_PHOTOMETRY_FIT_H
#define LSST_JOINTCAL_PHOTOMETRY_FIT_H

#include <string>
#include <iostream>
#include <sstream>

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/Tripletlist.h"
#include "lsst/jointcal/PhotometryModel.h"

namespace lsst {
namespace jointcal {

//! Class that handles the photometric least squares problem.
class PhotometryFit : public FitterBase {
private:
    // Associations &_associations;
    // std::string _whatToFit;
    bool _fittingModel, _fittingFluxes;
    PhotometryModel &_photometryModel;

    // counts in parameter subsets.
    unsigned int _nParModel;
    unsigned int _nParFluxes;
    // unsigned int _nParTot;

public:
    //! this is the only constructor
    PhotometryFit(Associations &associations, PhotometryModel &photometryModel)
            : FitterBase(associations), _photometryModel(photometryModel) {
        // The various _npar... are initialized in assignIndices.
        // Although there is no reason to adress them before one might be tempted by
        // evaluating a Chi2 rightaway, .. which uses these counts, so:
        assignIndices("");
    }

    /**
     * Set parameters to fit and assign indices in the big matrix.
     *
     * @param[in]  whatToFit  Valid strings : "Model", "Fluxes", which define
     *                        which parameter sets are going to be fitted.
     *                        whatToFit="Model Fluxes"  will set both parameter
     *                        sets variable when computing derivatives. Provided
     *                        it contains "Model", whatToFit is passed over to the
     *                        PhotometryModel, and can hence be used to control more
     *                        finely which subsets of the photometric model are
     *                        being fitted, if the the actual PhotometryModel
     *                        implements such a possibility.
     */
    void assignIndices(const std::string &whatToFit);

    /**
     * Offset the parameters by the requested quantities. The used parameter
     * layout is the one from the last call to assignIndices or minimize(). There
     * is no easy way to check that the current setting of whatToFit and the
     * provided Delta vector are compatible. We can only test the size.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    void offsetParams(const Eigen::VectorXd &delta);

    void saveResultTuples(const std::string &tupleName) const;

private:
    void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const;

    void accumulateStatRefStars(Chi2Accumulator &accum) const;

    void setMeasuredStarIndices(const MeasuredStar &measuredStar, std::vector<unsigned> &indices) const;

    /** Compute the derivatives of the measured stars and model for a CcdImage.
     *
     * The last argument allows to to process a sub-list for outlier removal.
     */
    void leastSquareDerivativesMeasurement(const CcdImage &ccdImage, TripletList &tripletList,
                                           Eigen::VectorXd &grad,
                                           const MeasuredStarList *measuredStarList = nullptr) const;

    /// Compute the derivatives of the reference terms
    void leastSquareDerivativesReference(const FittedStarList &fittedStarList, TripletList &tripletList,
                                         Eigen::VectorXd &grad) const;

#ifdef STORAGE
    //! Produces a tuple containing residuals of measurement terms.
    void makeMeasResTuple(const std::string &tupleName) const;

    //! Produces a tuple containing residuals of reference terms.
    void makeRefResTuple(const std::string &tupleName) const;

private:
    Point transformFittedStar(const FittedStar &fittedStar, const Gtransfo *sky2TP,
                              const Point &refractionVector, double refractionCoeff, double mjd) const;

    //! only for outlier removal
    void getMeasuredStarIndices(const MeasuredStar &measuredStar, std::vector<unsigned> &indices) const;
#endif
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_PHOTOMETRY_FIT_H
