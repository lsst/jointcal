#ifndef LSST_JOINTCAL_PHOTOMETRY_FIT_H
#define LSST_JOINTCAL_PHOTOMETRY_FIT_H

#include <string>
#include <iostream>
#include <sstream>

#include "lsst/log/Log.h"
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
public:
    //! this is the only constructor
    PhotometryFit(std::shared_ptr<Associations> associations,
                  std::shared_ptr<PhotometryModel> photometryModel)
            : FitterBase(associations),
              _fittingModel(false),
              _fittingFluxes(false),
              _photometryModel(photometryModel),
              _nParModel(0),
              _nParFluxes(0) {
        _log = LOG_GET("jointcal.PhotometryFit");
    }

    /// No copy or move: there is only ever one fitter of a given type.
    PhotometryFit(PhotometryFit const &) = delete;
    PhotometryFit(PhotometryFit &&) = delete;
    PhotometryFit &operator=(PhotometryFit const &) = delete;
    PhotometryFit &operator=(PhotometryFit &&) = delete;

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
    void assignIndices(std::string const &whatToFit) override;

    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc FitterBase::saveChi2MeasContributions
    void saveChi2MeasContributions(std::string const &baseName) const override;

    /// @copydoc FitterBase::saveChi2RefContributions
    void saveChi2RefContributions(std::string const &baseName) const override;

private:
    bool _fittingModel, _fittingFluxes;
    std::shared_ptr<PhotometryModel> _photometryModel;

    // counts in parameter subsets.
    unsigned int _nParModel;
    unsigned int _nParFluxes;

    void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const override;

    void accumulateStatRefStars(Chi2Accumulator &accum) const override;

    void getIndicesOfMeasuredStar(MeasuredStar const &measuredStar,
                                  std::vector<unsigned> &indices) const override;

    void leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                           Eigen::VectorXd &grad,
                                           MeasuredStarList const *measuredStarList = nullptr) const override;

    /// Compute the derivatives of the reference terms
    void leastSquareDerivativesReference(FittedStarList const &fittedStarList, TripletList &tripletList,
                                         Eigen::VectorXd &grad) const override;

#ifdef STORAGE
    Point transformFittedStar(FittedStar const &fittedStar, Gtransfo const *sky2TP,
                              Point const &refractionVector, double refractionCoeff, double mjd) const;
#endif
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_PHOTOMETRY_FIT_H
