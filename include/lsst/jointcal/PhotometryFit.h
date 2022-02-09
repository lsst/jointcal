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

#ifndef LSST_JOINTCAL_PHOTOMETRY_FIT_H
#define LSST_JOINTCAL_PHOTOMETRY_FIT_H

#include <iostream>
#include <sstream>
#include <string>
#include <utility>

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
    /**
     * Construct a photometry fitter.
     *
     * @param associations The associations catalog to use in the fitter.
     * @param photometryModel The model to build the fitter for.
     */
    PhotometryFit(std::shared_ptr<Associations> associations,
                  std::shared_ptr<PhotometryModel> photometryModel)
            : FitterBase(std::move(associations)),
              _fittingModel(false),
              _fittingFluxes(false),
              _photometryModel(std::move(photometryModel)) {
        _log = LOG_GET("lsst.jointcal.PhotometryFit");
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

    /// Return the model being fit.
    std::shared_ptr<PhotometryModel> getModel() const { return _photometryModel; }

protected:
    /// @copydoc FitterBase::saveChi2MeasContributions
    void saveChi2MeasContributions(std::string const &filename) const override;

    /// @copydoc FitterBase::saveChi2RefContributions
    void saveChi2RefContributions(std::string const &filename) const override;

private:
    bool _fittingModel, _fittingFluxes;
    std::shared_ptr<PhotometryModel> _photometryModel;

    void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const override;

    void accumulateStatRefStars(Chi2Accumulator &accum) const override;

    void getIndicesOfMeasuredStar(MeasuredStar const &measuredStar, IndexVector &indices) const override;

    void leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                           Eigen::VectorXd &grad,
                                           MeasuredStarList const *measuredStarList = nullptr) const override;

    /// Compute the derivatives of the reference terms
    void leastSquareDerivativesReference(FittedStarList const &fittedStarList, TripletList &tripletList,
                                         Eigen::VectorXd &grad) const override;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_PHOTOMETRY_FIT_H
