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

namespace lsst {
namespace jointcal {

void PhotometryFit::leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                                      Eigen::VectorXd &grad,
                                                      MeasuredStarList const *measuredStarList) const {
    /**********************************************************************/
    /* @note the math in this method and accumulateStatImageList() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* this routine works in two different ways: either providing the
       Ccd, of providing the MeasuredStarList. In the latter case, the
       Ccd should match the one(s) in the list. */
    if (measuredStarList) assert(&(measuredStarList->front()->getCcdImage()) == &ccdImage);

    unsigned nparModel = (_fittingModel) ? _photometryModel->getNpar(ccdImage) : 0;
    unsigned nparFlux = (_fittingFluxes) ? 1 : 0;
    unsigned nparTotal = nparModel + nparFlux;
    std::vector<unsigned> indices(nparModel, -1);
    if (_fittingModel) _photometryModel->getMappingIndices(ccdImage, indices);

    Eigen::VectorXd H(nparTotal);  // derivative matrix
    // current position in the Jacobian
    unsigned kTriplets = tripletList.getNextFreeIndex();
    const MeasuredStarList &catalog = (measuredStarList) ? *measuredStarList : ccdImage.getCatalogForFit();

    for (auto const &measuredStar : catalog) {
        if (!measuredStar->isValid()) continue;
// tweak the measurement errors
#ifdef FUTURE
        TweakPhotomMeasurementErrors(inPos, *measuredStar, _fluxError);
#endif
        H.setZero();  // we cannot be sure that all entries will be overwritten.

        double residual = _photometryModel->computeResidual(ccdImage, *measuredStar);
        double inverseSigma = 1.0 / _photometryModel->transformError(ccdImage, *measuredStar);
        double W = std::pow(inverseSigma, 2);

        if (_fittingModel) {
            _photometryModel->computeParameterDerivatives(*measuredStar, ccdImage, H);
            for (unsigned k = 0; k < indices.size(); k++) {
                unsigned l = indices[k];
                tripletList.addTriplet(l, kTriplets, H[k] * inverseSigma);
                grad[l] += H[k] * W * residual;
            }
        }
        if (_fittingFluxes) {
            unsigned index = measuredStar->getFittedStar()->getIndexInMatrix();
            // Note: H = dR/dFittedStarFlux == -1
            tripletList.addTriplet(index, kTriplets, -1.0 * inverseSigma);
            grad[index] += -1.0 * W * residual;
        }
        kTriplets += 1;  // each measurement contributes 1 column in the Jacobian
    }

    tripletList.setNextFreeIndex(kTriplets);
}

void PhotometryFit::leastSquareDerivativesReference(FittedStarList const &fittedStarList,
                                                    TripletList &tripletList, Eigen::VectorXd &grad) const {
    /**********************************************************************/
    /** @note the math in this method and accumulateStatReference() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    // Derivatives of terms involving fitted and refstars only contribute if we are fitting fluxes.
    if (!_fittingFluxes) return;
    // Can't compute anything if there are no refStars.
    if (_associations->refStarList.size() == 0) return;

    unsigned kTriplets = tripletList.getNextFreeIndex();

    for (auto const &fittedStar : fittedStarList) {
        auto refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;  // no contribution if no associated refstar
        // TODO: Can we actually work with multiple filters at a time in this scheme?
        // TODO: I feel we might only be able to fit one filter at a time, because these terms are
        // independent of the ccdImage, so we don't have a specific filter.
        // filter = ccdImage.getFilter();

        // W == inverseSigma^2
        double inverseSigma = 1.0 / refStar->getFluxErr();
        // Residual is fittedStar.flux - refStar.flux for consistency with measurement terms.
        double residual = fittedStar->getFlux() - refStar->getFlux();

        unsigned index = fittedStar->getIndexInMatrix();
        // Note: H = dR/dFittedStar == 1
        tripletList.addTriplet(index, kTriplets, 1.0 * inverseSigma);
        grad(index) += 1.0 * std::pow(inverseSigma, 2) * residual;
        kTriplets += 1;
    }
    tripletList.setNextFreeIndex(kTriplets);
}

void PhotometryFit::accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesMeasurement() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/
    for (auto const &ccdImage : ccdImageList) {
        auto &catalog = ccdImage->getCatalogForFit();

        for (auto const &measuredStar : catalog) {
            if (!measuredStar->isValid()) continue;
            double sigma = _photometryModel->transformError(*ccdImage, *measuredStar);
#ifdef FUTURE
            TweakPhotomMeasurementErrors(inPos, measuredStar, _fluxError);
#endif
            double residual = _photometryModel->computeResidual(*ccdImage, *measuredStar);

            double chi2Val = std::pow(residual / sigma, 2);
            accum.addEntry(chi2Val, 1, measuredStar);
        }  // end loop on measurements
    }
}

void PhotometryFit::accumulateStatRefStars(Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesReference() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    FittedStarList &fittedStarList = _associations->fittedStarList;
    for (auto const &fittedStar : fittedStarList) {
        auto refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;
        double chi2 = std::pow(((fittedStar->getFlux() - refStar->getFlux()) / refStar->getFluxErr()), 2);
        accum.addEntry(chi2, 1, fittedStar);
    }
}

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void PhotometryFit::getIndicesOfMeasuredStar(MeasuredStar const &measuredStar,
                                             std::vector<unsigned> &indices) const {
    indices.clear();
    if (_fittingModel) {
        _photometryModel->getMappingIndices(measuredStar.getCcdImage(), indices);
    }
    if (_fittingFluxes) {
        std::shared_ptr<FittedStar const> const fs = measuredStar.getFittedStar();
        unsigned fsIndex = fs->getIndexInMatrix();
        indices.push_back(fsIndex);
    }
}

void PhotometryFit::assignIndices(std::string const &whatToFit) {
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: now fitting: " << whatToFit);
    _fittingModel = (_whatToFit.find("Model") != std::string::npos);
    _fittingFluxes = (_whatToFit.find("Fluxes") != std::string::npos);
    // When entering here, we assume that whatToFit has already been interpreted.

    _nParModel = (_fittingModel) ? _photometryModel->assignIndices(whatToFit, 0) : 0;
    unsigned ipar = _nParModel;

    if (_fittingFluxes) {
        for (auto &fittedStar : _associations->fittedStarList) {
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in getIndicesOfMeasuredStar
            fittedStar->setIndexInMatrix(ipar);
            ipar += 1;
        }
    }
    _nParFluxes = ipar - _nParModel;
    _nParTot = ipar;
    LOGLS_DEBUG(_log,
                "nParameters total: " << _nParTot << " model: " << _nParModel << " fluxes: " << _nParFluxes);
}

void PhotometryFit::offsetParams(Eigen::VectorXd const &delta) {
    if (delta.size() != _nParTot)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "PhotometryFit::offsetParams : the provided vector length is not compatible with "
                          "the current whatToFit setting");
    if (_fittingModel) _photometryModel->offsetParams(delta);

    if (_fittingFluxes) {
        for (auto &fittedStar : _associations->fittedStarList) {
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when assigning indices (assignIndices())
            unsigned index = fittedStar->getIndexInMatrix();
            fittedStar->getFlux() -= delta(index);
        }
    }
}

void PhotometryFit::saveChi2MeasContributions(std::string const &baseName) const {
    std::ofstream ofile(baseName.c_str());
    std::string separator = "\t";

    ofile << "#id" << separator << "xccd" << separator << "yccd" << separator;
    ofile << "mag" << separator << "instMag" << separator << "instMagErr" << separator;
    ofile << "instFlux" << separator << "instFluxErr" << separator;
    ofile << "inputFlux" << separator << "inputFluxErr" << separator;
    ofile << "transformedFlux" << separator << "transformedFluxErr" << separator;
    ofile << "fittedFlux" << separator;
    ofile << "mjd" << separator << "color" << separator;
    ofile << "fsindex" << separator;
    ofile << "ra" << separator << "dec" << separator;
    ofile << "chi2" << separator << "nm" << separator;
    ofile << "chip" << separator << "visit" << separator << std::endl;

    ofile << "#id in source catalog" << separator << "coordinates in CCD" << separator << separator;
    ofile << "fitted magnitude" << separator;
    ofile << "measured instrument flux (ADU)" << separator << "measured instrument flux error" << separator;
    ofile << "measured flux (maggies)" << separator << "measured flux error" << separator;
    ofile << separator << separator;
    ofile << "fitted flux (maggies)" << separator;
    ofile << "modified Julian date of the measurement" << separator << "currently unused" << separator;
    ofile << "unique index of the fittedStar" << separator;
    ofile << "on-sky position of fitted star" << separator << separator;
    ofile << "contribution to Chi2 (1 dof)" << separator << "number of measurements of this FittedStar"
          << separator;
    ofile << "chip id" << separator << "visit id" << std::endl;

    const CcdImageList &ccdImageList = _associations->getCcdImageList();
    for (auto const &ccdImage : ccdImageList) {
        const MeasuredStarList &cat = ccdImage->getCatalogForFit();
        for (auto const &measuredStar : cat) {
            if (!measuredStar->isValid()) continue;
            double sigma = _photometryModel->transformError(*ccdImage, *measuredStar);
#ifdef FUTURE
            tweakPhotomMeasurementErrors(inPos, measuredStar, _fluxError);
#endif
            double flux = _photometryModel->transform(*ccdImage, *measuredStar);
            double fluxErr = _photometryModel->transformError(*ccdImage, *measuredStar);
            double jd = ccdImage->getMjd();
            std::shared_ptr<FittedStar const> const fittedStar = measuredStar->getFittedStar();
            double residual = _photometryModel->computeResidual(*ccdImage, *measuredStar);
            double chi2Val = std::pow(residual / sigma, 2);

            ofile << std::setprecision(9);
            ofile << measuredStar->getId() << separator << measuredStar->x << separator << measuredStar->y
                  << separator;
            ofile << fittedStar->getMag() << separator;
            ofile << measuredStar->getInstFlux() << separator << measuredStar->getInstFluxErr() << separator;
            ofile << measuredStar->getFlux() << separator << measuredStar->getFluxErr() << separator;
            ofile << flux << separator << fluxErr << separator << fittedStar->getFlux() << separator;
            ofile << jd << separator << fittedStar->color << separator;
            ofile << fittedStar->getIndexInMatrix() << separator;
            ofile << fittedStar->x << separator << fittedStar->y << separator;
            ofile << chi2Val << separator << fittedStar->getMeasurementCount() << separator;
            ofile << ccdImage->getCcdId() << separator << ccdImage->getVisit() << std::endl;
        }  // loop on measurements in image
    }      // loop on images
}

void PhotometryFit::saveChi2RefContributions(std::string const &baseName) const {
    std::ofstream ofile(baseName.c_str());
    std::string separator = "\t";

    ofile << "#ra" << separator << "dec " << separator;
    ofile << "mag" << separator << "color" << separator;
    ofile << "refFlux" << separator << "refFluxErr" << separator;
    ofile << "fittedFlux" << separator << "fittedFluxErr" << separator;
    ofile << "fsindex" << separator << "chi2" << separator << "nm" << std::endl;

    ofile << "#coordinates of fittedStar" << separator << separator;
    ofile << "magnitude" << separator << "currently unused" << separator;
    ofile << "refStar flux (maggies)" << separator << "refStar fluxErr" << separator;
    ofile << "fittedStar flux (maggies)" << separator << "fittedStar fluxErr" << separator;
    ofile << "unique index of the fittedStar" << separator << "refStar contribution to Chi2 (1 dof)"
          << separator << "number of measurements of this FittedStar" << std::endl;

    // The following loop is heavily inspired from PhotometryFit::computeChi2()
    const FittedStarList &fittedStarList = _associations->fittedStarList;
    for (auto const &fittedStar : fittedStarList) {
        const RefStar *refStar = fittedStar->getRefStar();
        if (refStar == nullptr) continue;

        double chi2 = std::pow(((fittedStar->getFlux() - refStar->getFlux()) / refStar->getFluxErr()), 2);

        ofile << std::setprecision(9);
        ofile << fittedStar->x << separator << fittedStar->y << separator;
        ofile << fittedStar->getMag() << separator << fittedStar->color << separator;
        ofile << refStar->getFlux() << separator << refStar->getFluxErr() << separator;
        ofile << fittedStar->getFlux() << separator << fittedStar->getFluxErr() << separator;
        ofile << fittedStar->getIndexInMatrix() << separator << chi2 << separator
              << fittedStar->getMeasurementCount() << std::endl;
    }  // loop on FittedStars
}

}  // namespace jointcal
}  // namespace lsst
