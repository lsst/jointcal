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

void PhotometryFit::LSDerivatives(TripletList &tripletList, Eigen::VectorXd &rhs) const {
    auto ccdImageList = _associations.getCcdImageList();
    for (auto const &ccdImage : ccdImageList) {
        LSDerivativesPerCcdImage(*ccdImage, tripletList, rhs);
    }
    LSDerivativesReference(_associations.fittedStarList, tripletList, rhs);
}

void PhotometryFit::LSDerivativesPerCcdImage(const CcdImage &ccdImage, TripletList &tripletList,
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

    Eigen::VectorXd h(npar_max);
    Eigen::VectorXd grad(npar_max);
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
        h.setZero();  // we cannot be sure that all entries will be overwritten.

        double pf = _photometryModel->photomFactor(ccdImage, measuredStar);
        auto fs = measuredStar.getFittedStar();

        double residual = measuredStar.getFlux() - pf * fs->getFlux();

        if (_fittingModel) {
            _photometryModel->setIndicesAndDerivatives(measuredStar, ccdImage, indices, h);
            for (unsigned k = 0; k < indices.size(); k++) {
                unsigned l = indices[k];
                tripletList.addTriplet(l, kTriplets, h[k] * fs->getFlux() / fluxErr);
                rhs[l] += h[k] * residual / sqr(fluxErr);
            }
        }
        if (_fittingFluxes) {
            unsigned index = fs->getIndexInMatrix();
            tripletList.addTriplet(index, kTriplets, pf / fluxErr);
            rhs[index] += residual * pf / sqr(fluxErr);
        }
        kTriplets += 1;  // each measurement contributes 1 column in the Jacobian
    }                    // end loop on measurements
    tripletList.setNextFreeIndex(kTriplets);
}

void PhotometryFit::LSDerivativesReference(const FittedStarList &fittedStarList, TripletList &tripletList,
                                           Eigen::VectorXd &rhs) const {
    // Derivatives of terms involving fitted and refstars only contribute if we are fitting fluxes.
    if (!_fittingFluxes) return;
    // Can't compute anything if there are no refStars.
    if (_associations.refStarList.size() == 0) return;

    unsigned kTriplets = tripletList.getNextFreeIndex();

    for (auto const &fittedStar : fittedStarList) {
        auto refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;

        kTriplets += 1;
    }
    tripletList.setNextFreeIndex(kTriplets);
}

// This is almost a selection of lines of LSDerivatives(CcdImage ...)
/* This routine is template because it is used
both with its first argument as "const CcdImageList &" and "CcdImageList &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ListType, class Accum>
void PhotometryFit::accumulateStat(ListType &listType, Accum &accum) const {
    for (auto &im : listType) {
        /**********************************************************************/
        /**  Changes in this routine should be reflected into LSDerivatives  */
        /**********************************************************************/
        auto &ccdIMage = *im;
        auto &catalog = ccdIMage.getCatalogForFit();

        for (auto const &measuredStar : catalog) {
            if (!measuredStar->isValid()) continue;
            // tweak the measurement errors
            double sigma = measuredStar->getFluxErr();
#ifdef FUTURE
            TweakPhotomMeasurementErrors(inPos, measuredStar, _fluxError);
#endif

            double pf = _photometryModel->photomFactor(ccdIMage, *measuredStar);
            auto fs = measuredStar->getFittedStar();
            double residual = measuredStar->getFlux() - pf * fs->getFlux();
            double chi2Val = sqr(residual / sigma);
            accum.addEntry(chi2Val, 1, measuredStar);
        }  // end loop on measurements
    }
}

//! for the list of images in the provided  association and the reference stars, if any
Chi2Statistic PhotometryFit::computeChi2() const {
    Chi2Statistic chi2;
    accumulateStat(_associations.getCcdImageList(), chi2);
    // so far, chi2.ndof contains the number of squares.
    // So, subtract here the number of parameters.
    chi2.ndof -= _nParTot;
    return chi2;
}

void PhotometryFit::outliersContributions(MeasuredStarList &outliers, TripletList &tripletList,
                                          Eigen::VectorXd &grad) {
    for (auto &outlier : outliers) {
        MeasuredStarList tmp;
        tmp.push_back(outlier);
        const CcdImage &ccdImage = outlier->getCcdImage();
        LSDerivativesPerCcdImage(ccdImage, tripletList, grad, &tmp);
        outlier->setValid(false);
        auto fs = std::const_pointer_cast<FittedStar>(outlier->getFittedStar());
        fs->getMeasurementCount()--;
    }
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

void PhotometryFit::findOutliers(double nSigmaCut, MeasuredStarList &outliers) const {
    /* Aims at providing an outlier list for small-rank update
       of the factorization. */

    // collect chi2 contributions of the measurements.
    Chi2List chi2List;
    // chi2List.reserve(_nMeasuredStars);
    accumulateStat(_associations.ccdImageList, chi2List);
    // do some stat
    unsigned nval = chi2List.size();
    if (nval == 0) return;
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
            outliers.push_back(std::dynamic_pointer_cast<MeasuredStar>(chi2->star));
        }
    }  // end loop on measurements
    LOGLS_INFO(_log, "findMeasOutliers: found " << outliers.size() << " outliers");
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
    LSDerivatives(tripletList, grad);
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
                              << hessian.rows() << " nnz=" << hessian.nonZeros()
                              << " filling-frac = " << hessian.nonZeros() / sqr(hessian.rows()));

    Eigen::SimplicialLDLT<SpMat> chol(hessian);
    // CholmodSimplicialLDLT2<SpMat> chol(hessian);

    if (chol.info() != Eigen::Success) {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return 1;
    }

    Eigen::VectorXd delta = chol.solve(grad);

    offsetParams(delta);
    return returnCode;
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
            double pf = _photometryModel->photomFactor(im, ms);
            double jd = im.getMjd();
            auto fs = ms.getFittedStar();
            double residual = ms.getFlux() - pf * fs->getFlux();
            double chi2Val = sqr(residual / sigma);
            tuple << ms.x << ' ' << ms.y << ' ' << fs->getMag() << ' ' << ms.getFlux() << ' '
                  << ms.getFluxErr() << ' ' << fs->getFlux() << ' ' << pf << ' ' << jd << ' ' << fs->color
                  << ' ' << fs->getIndexInMatrix() << ' ' << fs->x << ' ' << fs->y << ' ' << chi2Val << ' '
                  << fs->getMeasurementCount() << ' ' << im.getCcdId() << ' ' << im.getVisit() << endl;
        }  // loop on measurements in image
    }      // loop on images
}
}  // namespace jointcal
}  // namespace lsst
