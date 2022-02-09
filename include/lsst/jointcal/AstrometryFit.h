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

#ifndef LSST_JOINTCAL_ASTROMETRY_FIT_H
#define LSST_JOINTCAL_ASTROMETRY_FIT_H

#include <string>
#include <iostream>
#include <sstream>

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/Tripletlist.h"
#include "lsst/jointcal/AstrometryModel.h"

namespace lsst {
namespace jointcal {

/**
 * Class that handles the astrometric least squares problem.
 *
 * This is the class that actually computes the quantities required to carry out
 * a LS astrometric fit wrt distortion mappings and coordinates of common
 * objects. Namely it computes the Jacobian and gradient of the chi2 (w.r.t.
 * parameters), and the Chi2 itself. It interfaces with the actual modelling of
 * distortions via a mimimum virtual interface AstrometryModel, and the actual
 * mappings via an other virtual interface : Mapping.
 *
 * In short AstrometryFit aims at computing derivatives of least quares. The terms
 * of the chi2 are of two kinds:
 *
 * kind 1 ->   (T(X_M) - p(F))^T W (T(X_M) - p(F))
 *
 * with X_M is a measured (2d) position in CCD coordinates, F refers to the
 * position of the object in some space, defined in practise by p. There is one
 * such term per measurement. The default setup would be that p is the
 * projection from sky to some tangent plane and hence T maps the CCD
 * coordinates onto this TP. p is obtained via the DistorsionModel and can be
 * different for all CcdImage's. Depending on what is beeing fitted, one could
 * imagine cases where the projector p is the same for all CcdImages.
 *
 * Kind 2  -> (p'(F)-p'(R))^T  W_R (p'(F)-p'(R)) R refers to some
 * externally-provided reference object position, and p' to some projector from
 * sky to some plane. The reference objects define the overall coordinate frame,
 * which is required when all T and all F are fitted  simultaneously. There is
 * one such term per external reference object. There can be more F (fitted)
 * objects than R (reference) objects.
 *
 * In the same framework, one can fit relative transforms between images by
 * setting p = Identity for all input CcdImages and not fitting T for one of the
 * CcdImage's. One does not need reference object and would then naturally not
 * have any Kind 2 terms.
 */
class AstrometryFit : public FitterBase {
public:
    //! this is the only constructor
    AstrometryFit(std::shared_ptr<Associations> associations,
                  std::shared_ptr<AstrometryModel> astrometryModel, double posError);

    /// No copy or move: there is only ever one fitter of a given type.
    AstrometryFit(AstrometryFit const &) = delete;
    AstrometryFit(AstrometryFit &&) = delete;
    AstrometryFit &operator=(AstrometryFit const &) = delete;
    AstrometryFit &operator=(AstrometryFit &&) = delete;

    /**
     * Set parameters to fit and assign indices in the big matrix.
     *
     * @param[in]  whatToFit   Valid strings: zero or more of "Distortions", "Positions",
     *                         "PM" which define which parameter set
     *                         is going to be variable when computing
     *                         derivatives (leastSquareDerivatives) and minimizing
     *                         (minimize()). whatToFit="Positions Distortions"
     *                         will minimize w.r.t mappings and objects
     *                         positions, and not w.r.t proper motions and
     *                         refraction modeling. However if proper motions
     *                         and/or refraction parameters have already been
     *                         set, then they are accounted for when computing
     *                         residuals.  The string is forwarded to the
     *                         AstrometryModel, and it can then be used to turn
     *                         subsets of distortion parameter on or off, if the
     *                         AstrometryModel implements such a thing.
     */
    void assignIndices(std::string const &whatToFit) override;

    /**
     * The transformations used to propagate errors are freezed to the current
     * state. The routine can be called when the mappings are roughly in place.
     * After the call, the transformations used to propage errors are no longer
     * affected when updating the mappings. This allows to have an exactly linear
     * fit, which can be useful.
     */
    void freezeErrorTransform() { _astrometryModel->freezeErrorTransform(); }

    void offsetParams(Eigen::VectorXd const &delta) override;

    /// Return the model being fit.
    std::shared_ptr<AstrometryModel> getModel() const { return _astrometryModel; }

    /**
     * DEBUGGING routine
     */
    void checkStuff();

protected:
    /// @copydoc FitterBase::saveChi2MeasContributions
    void saveChi2MeasContributions(std::string const &filename) const override;

    /// @copydoc FitterBase::saveChi2RefContributions
    void saveChi2RefContributions(std::string const &filename) const override;

private:
    bool _fittingDistortions{}, _fittingPos{}, _fittingPM{};
    std::shared_ptr<AstrometryModel> _astrometryModel;

    double _epoch;     // epoch to correct proper motion/parallax to (Julian Epoch year, e.g. J2000.0)
    double _posError;  // constant term on error on position (in pixel unit)

    void leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                           Eigen::VectorXd &grad,
                                           MeasuredStarList const *msList = nullptr) const override;

    void leastSquareDerivativesReference(FittedStarList const &fittedStarList, TripletList &tripletList,
                                         Eigen::VectorXd &grad) const override;

    void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const override;

    void accumulateStatRefStars(Chi2Accumulator &accum) const override;

    void getIndicesOfMeasuredStar(MeasuredStar const &measuredStar, IndexVector &indices) const override;

    /**
     * Transform the positions of a FittedStar into the frame of a MeasuredStar.
     *
     * @param fittedStar The star to transform.
     * @param sky2TP Transformation from sky coordinates to CcdImage tangent plane.
     * @param deltaYears Difference in years between FittedStar and MeasuredStar epochs.
     * @return Corrected position of FittedStar.
     */
    Point transformFittedStar(FittedStar const &fittedStar, AstrometryTransform const &sky2TP,
                              double deltaYears) const;

    /// Compute the chi2 (per star or total, depending on which Chi2Accumulator is used) from one CcdImage.
    void accumulateStatImage(CcdImage const &ccdImage, Chi2Accumulator &accum) const;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_ASTROMETRY_FIT_H
