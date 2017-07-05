#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include "Eigen/Sparse"

#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/PhotometryFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Tripletlist.h"

using namespace std;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.PhotometryFit");
}

namespace lsst {
namespace jointcal {

PhotometryFit::PhotometryFit(Associations &associations, PhotometryModel *photometryModel)
        : _associations(associations), _photometryModel(photometryModel), _lastNTrip(0) {
    // The various _npar... are initialized in assignIndices.
    // Although there is no reason to adress them before one might be tempted by
    // evaluating a Chi2 rightaway, .. which uses these counts, so:
    assignIndices("");
}

void PhotometryFit::leastSquareDerivatives(TripletList &tripletList, Eigen::VectorXd &rhs) const {
    auto ccdImageList = _associations.getCcdImageList();
    for (auto const &ccdImage : ccdImageList) {
        leastSquareDerivativesMeasurement(*ccdImage, tripletList, rhs);
    }
    leastSquareDerivativesReference(_associations.fittedStarList, tripletList, rhs);
}

void PhotometryFit::leastSquareDerivativesMeasurement(const CcdImage &ccdImage, TripletList &tripletList,
                                                      Eigen::VectorXd &rhs,
                                                      const MeasuredStarList *measuredStarList) const {
    /***************************************************************************/
    /**  Changes in this routine should be reflected into accumulateStat       */
    /***************************************************************************/
    /* this routine works in two different ways: either providing the
       Ccd, of providing the MeasuredStarList. In the latter case, the
       Ccd should match the one(s) in the list. */
    if (measuredStarList) assert(&(measuredStarList->front()->getCcdImage()) == &ccdImage);

    unsigned npar_max = 100;  // anything large
    vector<unsigned> indices(npar_max, -1);

    Eigen::VectorXd H(npar_max);
    // current position in the Jacobian
    unsigned kTriplets = tripletList.getNextFreeIndex();
    const MeasuredStarList &catalog = (measuredStarList) ? *measuredStarList : ccdImage.getCatalogForFit();

    for (auto const &i : catalog) {
        const MeasuredStar &measuredStar = *i;
        if (!measuredStar.isValid()) continue;
        // tweak the measurement errors
        double fluxErr = measuredStar.getFluxErr();
#ifdef FUTURE
        TweakPhotomMeasurementErrors(inPos, measuredStar, _fluxError);
#endif
        H.setZero();  // we cannot be sure that all entries will be overwritten.

        double photomFactor = _photometryModel->photomFactor(ccdImage, measuredStar);
        auto fs = measuredStar.getFittedStar();

        double residual = measuredStar.getFlux() - photomFactor * fs->getFlux();

        if (_fittingModel) {
            _photometryModel->setIndicesAndDerivatives(measuredStar, ccdImage, indices, H);
            for (unsigned k = 0; k < indices.size(); k++) {
                unsigned l = indices[k];
                tripletList.addTriplet(l, kTriplets, H[k] * fs->getFlux() / fluxErr);
                rhs[l] += H[k] * residual / sqr(fluxErr);
            }
        }
        if (_fittingFluxes) {
            unsigned index = fs->getIndexInMatrix();
            tripletList.addTriplet(index, kTriplets, photomFactor / fluxErr);
            rhs[index] += residual * photomFactor / sqr(fluxErr);
        }
        kTriplets += 1;  // each measurement contributes 1 column in the Jacobian
    }                    // end loop on measurements
    tripletList.setNextFreeIndex(kTriplets);
}

void PhotometryFit::leastSquareDerivativesReference(const FittedStarList &fittedStarList,
                                                    TripletList &tripletList, Eigen::VectorXd &rhs) const {
    // Derivatives of terms involving fitted and refstars only contribute if we are fitting fluxes.
    if (!_fittingFluxes) return;
    // Can't compute anything if there are no refStars.
    if (_associations.refStarList.size() == 0) return;

    unsigned kTriplets = tripletList.getNextFreeIndex();

    for (auto const &fittedStar : fittedStarList) {
        auto refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;  // no contribution if no associated refstar
        // TODO: Can we actually work with multiple filters at a time in this scheme?
        // TODO: I feel we might only be able to fit one filter at a time, because these terms are
        // independent of the ccdImage, so we don't have a specific filter.
        // filter = ccdImage.getFilter();

        // the derivative of the residual of one refStar term
        double h = refStar->getFlux() - fittedStar->getFlux();
        // the flux weight of one refStar term
        double w = 1.0 / refStar->getFluxErr();
        // Residual is fittedStar.flux - refStar.flux for consistency with measurement terms.
        double residual = fittedStar->getFlux() - refStar->getFlux();

        // Because we have no projector here (unlike AstrometryFit), the calculation of the
        // error contribution due to the reference stars is much simpler.
        unsigned index = fittedStar->getIndexInMatrix();
        tripletList.addTriplet(index, kTriplets, h * w);
        rhs(index) += residual * sqr(refStar->getFluxErr());
        kTriplets += 1;
    }
    tripletList.setNextFreeIndex(kTriplets);
}

// This is almost a selection of lines of leastSquareDerivatives(CcdImage ...)
/* This routine is template because it is used
both with its first argument as "const CcdImageList &" and "CcdImageList &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ListType, class Accum>
void PhotometryFit::accumulateStatImageList(ListType &listType, Accum &accum) const {
    for (auto ccdImage : listType) {
        /**********************************************************************/
        /**  Changes in this routine should be reflected into leastSquareDerivatives  */
        /**********************************************************************/
        auto &catalog = ccdImage->getCatalogForFit();

        for (auto const &measuredStar : catalog) {
            if (!measuredStar->isValid()) continue;
            // tweak the measurement errors
            double sigma = measuredStar->getFluxErr();
#ifdef FUTURE
            TweakPhotomMeasurementErrors(inPos, measuredStar, _fluxError);
#endif

            double photomFactor = _photometryModel->photomFactor(*ccdImage, *measuredStar);
            auto fs = measuredStar->getFittedStar();
            double residual = measuredStar->getFlux() - photomFactor * fs->getFlux();
            double chi2Val = sqr(residual / sigma);
            accum.addEntry(chi2Val, 1, measuredStar);
        }  // end loop on measurements
    }
}

template <class Accum>
void PhotometryFit::accumulateStatRefStars(Accum &accum) const {
    FittedStarList &fittedStarList = _associations.fittedStarList;
    for (auto const &fittedStar : fittedStarList) {
        auto refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;
        double chi2 = sqr((fittedStar->getFlux() - refStar->getFlux()) / refStar->getFluxErr());
        accum.addEntry(chi2, 1, fittedStar);
    }
}

//! for the list of images in the provided  association and the reference stars, if any
Chi2Statistic PhotometryFit::computeChi2() const {
    Chi2Statistic chi2;
    accumulateStatImageList(_associations.getCcdImageList(), chi2);
    accumulateStatRefStars(chi2);
    // so far, chi2.ndof contains the number of squares.
    // So, subtract here the number of parameters.
    chi2.ndof -= _nParTot;
    return chi2;
}

void PhotometryFit::outliersContributions(MeasuredStarList &msOutliers, FittedStarList &fsOutliers,
                                          TripletList &tripletList, Eigen::VectorXd &grad) {
    for (auto &outlier : msOutliers) {
        MeasuredStarList tmp;
        tmp.push_back(outlier);
        const CcdImage &ccdImage = outlier->getCcdImage();
        leastSquareDerivativesMeasurement(ccdImage, tripletList, grad, &tmp);
    }
    leastSquareDerivativesReference(fsOutliers, tripletList, grad);
}

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void PhotometryFit::setMeasuredStarIndices(const MeasuredStar &measuredStar,
                                           std::vector<unsigned> &indices) const {
    indices.clear();
    if (_fittingModel) {
        Eigen::VectorXd h(100);
        _photometryModel->setIndicesAndDerivatives(measuredStar, measuredStar.getCcdImage(), indices, h);
    }
    if (_fittingFluxes) {
        auto fs = measuredStar.getFittedStar();
        unsigned fsIndex = fs->getIndexInMatrix();
        indices.push_back(fsIndex);
    }
}

unsigned PhotometryFit::findOutliers(double nSigmaCut, MeasuredStarList &msOutliers,
                                     FittedStarList &fsOutliers) const {
    /* Aims at providing an outlier list for small-rank update
       of the factorization. */

    // collect chi2 contributions of the measurements.
    Chi2List chi2List;
    // chi2List.reserve(_nMeasuredStars);
    accumulateStatImageList(_associations.ccdImageList, chi2List);
    accumulateStatRefStars(chi2List);
    // do some stat
    unsigned nval = chi2List.size();
    if (nval == 0) return 0;
    sort(chi2List.begin(), chi2List.end());  // increasing order. We rely on this later.
    double median = (nval & 1) ? chi2List[nval / 2].chi2
                               : 0.5 * (chi2List[nval / 2 - 1].chi2 + chi2List[nval / 2].chi2);
    // some more stats. should go into the class if recycled anywhere else
    auto averageAndSigma = chi2List.computeAverageAndSigma();
    LOGLS_DEBUG(_log, "RemoveOutliers chi2 stat: mean/median/sigma " << averageAndSigma.first << '/' << median
                                                                     << '/' << averageAndSigma.second);
    double cut = averageAndSigma.first + nSigmaCut * averageAndSigma.second;
    /* For each of the parameters, we will not remove more than 1
       measurement that contributes to constraining it. Keep track
       of the affected parameters using an integer vector. This is the
       trick that Marc Betoule came up to for outlier removals in "star
       flats" fits. */
    Eigen::VectorXi affectedParams(_nParTot);
    affectedParams.setZero();

    unsigned nOutliers = 0;
    // start from the strongest outliers, i.e. at the end of the array.
    for (auto chi2 = chi2List.rbegin(); chi2 != chi2List.rend(); ++chi2) {
        if (chi2->chi2 < cut) break;  // because the array is sorted.
        vector<unsigned> indices;
        setMeasuredStarIndices(*(chi2->star), indices);
        bool drop_it = true;
        /* find out if a stronger outlier constraining one of the parameters
           this one contrains was already discarded. If yes, we keep this one */
        for (auto const &chi2 : indices)
            if (affectedParams(chi2) != 0) drop_it = false;

        if (drop_it) {
            for (auto const &i : indices) affectedParams(i)++;
            msOutliers.push_back(std::dynamic_pointer_cast<MeasuredStar>(chi2->star));
        }
        nOutliers++;
    }  // end loop on measurements
    LOGLS_INFO(_log, "findOutliers: found " << msOutliers.size() << " meas outliers and " << fsOutliers.size()
                                            << " ref outliers ");
    return nOutliers;
}

void PhotometryFit::assignIndices(const std::string &whatToFit) {
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: now fitting: " << whatToFit);
    _fittingModel = (_whatToFit.find("Model") != string::npos);
    _fittingFluxes = (_whatToFit.find("Fluxes") != string::npos);
    // When entering here, we assume that whatToFit has already been interpreted.

    _nParModel = (_fittingModel) ? _photometryModel->assignIndices(whatToFit, 0) : 0;
    unsigned ipar = _nParModel;

    if (_fittingFluxes) {
        FittedStarList &fsl = _associations.fittedStarList;
        for (auto &i : fsl) {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in SetMeasuredStarIndices
            fs.setIndexInMatrix(ipar);
            ipar += 1;
        }
    }
    _nParFluxes = ipar - _nParModel;
    _nParTot = ipar;
}

void PhotometryFit::offsetParams(const Eigen::VectorXd &delta) {
    if (delta.size() != _nParTot)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "PhotometryFit::offsetParams : the provided vector length is not compatible with "
                          "the current whatToFit setting");
    if (_fittingModel) _photometryModel->offsetParams(delta);

    if (_fittingFluxes) {
        FittedStarList &fsl = _associations.fittedStarList;
        for (auto &i : fsl) {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when assigning indices (assignIndices())
            unsigned index = fs.getIndexInMatrix();
            fs.getFlux() += delta(index);
        }
    }
}

int PhotometryFit::minimize(const std::string &whatToFit, double nSigmaCut) {
    assignIndices(whatToFit);

    // return code can take 3 values :
    // 0 : fit has converged - no more outliers
    // 1 : still some ouliers but chi2 increases
    // 2 : factorization failed
    int returnCode = 0;

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
                              << " filling-frac = " << hessian.nonZeros() / sqr(hessian.rows()));

    Eigen::SimplicialLDLT<SpMat> chol(hessian);
    // CholmodSimplicialLDLT2<SpMat> chol(hessian);

    if (chol.info() != Eigen::Success) {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return 2;
    }

    double old_chi2 = computeChi2().chi2;

    Eigen::VectorXd delta = chol.solve(grad);
    offsetParams(delta);
    Chi2Statistic current_chi2(computeChi2());
    LOGLS_DEBUG(_log, current_chi2);
    if (current_chi2.chi2 > old_chi2) {
        LOGL_WARN(_log, "chi2 went up, skipping outlier rejection");
        returnCode = 1;
        return returnCode;
    }

    if (nSigmaCut == 0) return returnCode;

    MeasuredStarList msOutliers;
    FittedStarList fsOutliers;
    int n_outliers = findOutliers(nSigmaCut, msOutliers, fsOutliers);
    if (n_outliers == 0) return returnCode;
    outliersContributions(msOutliers, fsOutliers, tripletList, grad);
    removeMeasOutliers(msOutliers);
    removeRefOutliers(fsOutliers);
    LOGLS_INFO(_log, "Total number of outliers " << n_outliers);

    return returnCode;
}

void PhotometryFit::removeMeasOutliers(MeasuredStarList &outliers) {
    for (auto &measuredStar : outliers) {
        auto fittedStar = std::const_pointer_cast<FittedStar>(measuredStar->getFittedStar());
        measuredStar->setValid(false);
        fittedStar->getMeasurementCount()--;  // could be put in setValid
    }
}

void PhotometryFit::removeRefOutliers(FittedStarList &outliers) {
    for (auto &fittedStar : outliers) {
        fittedStar->setRefStar(nullptr);
    }
}

void PhotometryFit::makeResTuple(const std::string &tupleName) const {
    std::ofstream tuple(tupleName.c_str());
    /* If we think the some coordinate on the focal plane is relevant in
       the ntuple, because thmodel relies on it, then we have to add
       some function to the model that returns this relevant
       coordinate. */
    tuple << "#xccd: coordinate in CCD" << endl
          << "#yccd: " << endl
          << "#mag: rough mag" << endl
          << "#flux : measured flux" << endl
          << "#fluxError : measured flux error" << endl
          << "#fflux : fitted flux" << endl
          << "#phot_factor:" << endl
          << "#jd: Julian date of the measurement" << endl
          << "#color : " << endl
          << "#fsindex: some unique index of the object" << endl
          << "#ra: pos of fitted star" << endl
          << "#dec: pos of fitted star" << endl
          << "#chi2: contribution to Chi2 (1 dof)" << endl
          << "#nm: number of measurements of this FittedStar" << endl
          << "#chip: chip number" << endl
          << "#visit: visit id" << endl
          << "#end" << endl;
    const CcdImageList &ccdImageList = _associations.getCcdImageList();
    for (auto const &i : ccdImageList) {
        const CcdImage &im = *i;
        const MeasuredStarList &cat = im.getCatalogForFit();
        for (auto const &is : cat) {
            const MeasuredStar &ms = *is;
            if (!ms.isValid()) continue;
            double sigma = ms.getFluxErr();
#ifdef FUTURE
            tweakPhotomMeasurementErrors(inPos, ms, _fluxError);
#endif
            double photomFactor = _photometryModel->photomFactor(im, ms);
            double jd = im.getMjd();
            auto fs = ms.getFittedStar();
            double residual = ms.getFlux() - photomFactor * fs->getFlux();
            double chi2Val = sqr(residual / sigma);
            tuple << ms.x << ' ' << ms.y << ' ' << fs->getMag() << ' ' << ms.getFlux() << ' '
                  << ms.getFluxErr() << ' ' << fs->getFlux() << ' ' << photomFactor << ' ' << jd << ' '
                  << fs->color << ' ' << fs->getIndexInMatrix() << ' ' << fs->x << ' ' << fs->y << ' '
                  << chi2Val << ' ' << fs->getMeasurementCount() << ' ' << im.getCcdId() << ' '
                  << im.getVisit() << endl;
        }  // loop on measurements in image
    }      // loop on images
}
}  // namespace jointcal
}  // namespace lsst
