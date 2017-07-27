// -*- LSST-C++ -*-
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
private:
    bool _fittingDistortions, _fittingPos, _fittingRefrac, _fittingPM;
    AstrometryModel &_astrometryModel;
    double _referenceColor, _sigCol;  // average and r.m.s color
    unsigned _nRefrac;
    double _refractionCoefficient;    // fit parameter
    unsigned int _refracPosInMatrix;  // where it stands
    double _JDRef;                    // average Julian date

    // counts in parameter subsets.
    unsigned int _nParDistortions;
    unsigned int _nParPositions;
    double _posError;  // constant term on error on position (in pixel unit)

public:
    //! this is the only constructor
    AstrometryFit(Associations &associations, AstrometryModel &astrometryModel, double posError);

    /**
     * Set parameters to fit and assign indices in the big matrix.
     *
     * @param[in]  whatToFit   Valid strings: zero or more of "Distortions", "Positions",
     *                         "Refrac", "PM" which define which parameter set
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
    void assignIndices(const std::string &whatToFit);

    /**
     * The transformations used to propagate errors are freezed to the current
     * state. The routine can be called when the mappings are roughly in place.
     * After the call, the transformations used to propage errors are no longer
     * affected when updating the mappings. This allows to have an exactly linear
     * fit, which can be useful.
     */
    void freezeErrorScales() { _astrometryModel.freezeErrorScales(); }

    /**
     * Offset the parameters by the requested quantities. The used parameter
     * layout is the one from the last call to assignIndices or minimize().
     * There is no easy way to check that the current setting of whatToFit and
     * the provided Delta vector are compatible. We can only test the size.
     *
     * @param[in]  delta  vector of offsets to apply
     */
    void offsetParams(const Eigen::VectorXd &delta);

    void saveResultTuples(const std::string &tupleName) const;

    //! Produces a tuple containing residuals of measurement terms.
    void makeMeasResTuple(const std::string &tupleName) const;

    //! Produces a tuple containing residuals of reference terms.
    void makeRefResTuple(const std::string &tupleName) const;

    /**
     * DEBUGGING routine
     */
    void checkStuff();

private:
    // void leastSquareDerivatives(TripletList &tripletList, Eigen::VectorXd &grad) const;

    void leastSquareDerivativesMeasurement(const CcdImage &ccdImage, TripletList &tripletList,
                                           Eigen::VectorXd &grad,
                                           const MeasuredStarList *msList = nullptr) const;

    void leastSquareDerivativesReference(const FittedStarList &fittedStarList, TripletList &tripletList,
                                         Eigen::VectorXd &grad) const;

    Point transformFittedStar(const FittedStar &fittedStar, const Gtransfo *sky2TP,
                              const Point &refractionVector, double refractionCoeff, double mjd) const;

    void accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const;

    void accumulateStatImage(CcdImage const &ccdImage, Chi2Accumulator &accum) const;

    void accumulateStatRefStars(Chi2Accumulator &accum) const;

    //! only for outlier removal
    void setMeasuredStarIndices(const MeasuredStar &ms, std::vector<unsigned> &indices) const;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_ASTROMETRY_FIT_H
