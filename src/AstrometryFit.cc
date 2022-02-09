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

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>

#include "Eigen/Sparse"

#include "lsst/log/Log.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/AstrometryFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Chi2.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Tripletlist.h"

namespace {
/**
 * Compute the Chi2 of a refstar projected onto the fittedStar tangent point.
 *
 * The x/y position deltas are technically (refStar.x - fittedStar.x), but the latter is 0 because of the
 * tangent plane definition.
 *
 * @param refStar star projected onto the fittedStar tangent point.
 * @return chi2 contribution from this star.
 */
double computeProjectedRefStarChi2(lsst::jointcal::FatPoint refStar) {
    double det = refStar.vx * refStar.vy - std::pow(refStar.vxy, 2);
    double wxx = refStar.vy / det;
    double wyy = refStar.vx / det;
    double wxy = -refStar.vxy / det;
    return wxx * std::pow(refStar.x, 2) + 2 * wxy * refStar.x * refStar.y + wyy * std::pow(refStar.y, 2);
}
}  // namespace

namespace lsst {
namespace jointcal {

AstrometryFit::AstrometryFit(std::shared_ptr<Associations> associations,
                             std::shared_ptr<AstrometryModel> astrometryModel, double posError)
        : FitterBase(associations),
          _astrometryModel(std::move(astrometryModel)),
          _epoch(_associations->getEpoch()),
          _posError(posError) {
    _log = LOG_GET("lsst.jointcal.AstrometryFit");
}

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstrometryFit::transformFittedStar(FittedStar const &fittedStar, AstrometryTransform const &sky2TP,
                                         double deltaYears) const {
    Point fittedStarInTP;
    if (fittedStar.getRefStar()) {
        // PM data is on-sky, not in the tangent plane
        Point temp = fittedStar.getRefStar()->applyProperMotion(fittedStar, deltaYears);
        fittedStarInTP = sky2TP.apply(temp);
    } else {
        fittedStarInTP = sky2TP.apply(fittedStar);
    }
    return fittedStarInTP;
}

/*! This is the first implementation of an error "model".  We'll
  certainly have to upgrade it. MeasuredStar provides the mag in case
  we need it.  */
static void tweakAstromMeasurementErrors(FatPoint &P, MeasuredStar const &Ms, double error) {
    static bool called = false;
    static double increment = 0;
    if (!called) {
        increment = std::pow(error, 2);  // was in Preferences
        called = true;
    }
    P.vx += increment;
    P.vy += increment;
}

// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstrometryFit::leastSquareDerivativesMeasurement(CcdImage const &ccdImage, TripletList &tripletList,
                                                      Eigen::VectorXd &fullGrad,
                                                      MeasuredStarList const *msList) const {
    /**********************************************************************/
    /* @note the math in this method and accumulateStatImage() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* Setup */
    /* this routine works in two different ways: either providing the
       ccdImage, of providing the MeasuredStarList. In the latter case, the
       ccdImage should match the one(s) in the list. */
    if (msList) assert(&(msList->front()->getCcdImage()) == &ccdImage);

    // get the Mapping
    const AstrometryMapping *mapping = _astrometryModel->getMapping(ccdImage);
    // count parameters
    std::size_t npar_mapping = (_fittingDistortions) ? mapping->getNpar() : 0;
    std::size_t npar_pos = (_fittingPos) ? 2 : 0;
    std::size_t npar_pm = (_fittingPM) ? 2 : 0;
    std::size_t npar_tot = npar_mapping + npar_pos + npar_pm;
    // if (npar_tot == 0) this CcdImage does not contribute
    // any constraint to the fit, so :
    if (npar_tot == 0) return;
    IndexVector indices(npar_tot, -1);
    if (_fittingDistortions) mapping->getMappingIndices(indices);

    // FittedStar is "observed" epoch, MeasuredStar is "baseline"
    double deltaYears = _epoch - ccdImage.getEpoch();
    // transformation from sky to TP
    auto sky2TP = _astrometryModel->getSkyToTangentPlane(ccdImage);
    // reserve matrices once for all measurements
    AstrometryTransformLinear dypdy;
    // the shape of H (et al) is required this way in order to be able to
    // separate derivatives along x and y as vectors.
    Eigen::MatrixX2d H(npar_tot, 2), halpha(npar_tot, 2), HW(npar_tot, 2);
    Eigen::Matrix2d transW(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::VectorXd grad(npar_tot);
    // current position in the Jacobian
    Eigen::Index kTriplets = tripletList.getNextFreeIndex();
    const MeasuredStarList &catalog = (msList) ? *msList : ccdImage.getCatalogForFit();

    for (auto &i : catalog) {
        const MeasuredStar &ms = *i;
        if (!ms.isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = ms;
        tweakAstromMeasurementErrors(inPos, ms, _posError);
        H.setZero();  // we cannot be sure that all entries will be overwritten.
        FatPoint outPos;
        // should *not* fill H if whatToFit excludes mapping parameters.
        if (_fittingDistortions)
            mapping->computeTransformAndDerivatives(inPos, outPos, H);
        else
            mapping->transformPosAndErrors(inPos, outPos);

        std::size_t ipar = npar_mapping;
        double det = outPos.vx * outPos.vy - std::pow(outPos.vxy, 2);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, "Inconsistent measurement errors: drop measurement at "
                                     << Point(ms) << " in image " << ccdImage.getName());
            continue;
        }
        transW(0, 0) = outPos.vy / det;
        transW(1, 1) = outPos.vx / det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy / det;
        // compute alpha, a triangular square root
        // of transW (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(transW(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = transW(0, 1) / alpha(0, 0);
        // DB - I think that the next line is equivalent to : alpha(1,1) = 1./sqrt(outPos.vy)
        // PA - seems correct !
        alpha(1, 1) = 1. / sqrt(det * transW(0, 0));
        alpha(0, 1) = 0;

        std::shared_ptr<FittedStar const> const fs = ms.getFittedStar();

        Point fittedStarInTP = transformFittedStar(*fs, *sky2TP, deltaYears);

        // compute derivative of TP position w.r.t sky position ....
        if (npar_pos > 0)  // ... if actually fitting FittedStar position
        {
            sky2TP->computeDerivative(*fs, dypdy, 1e-3);
            // sign checked
            // TODO Still have to check with non trivial non-diagonal terms
            H(npar_mapping, 0) = -dypdy.A11();
            H(npar_mapping + 1, 0) = -dypdy.A12();
            H(npar_mapping, 1) = -dypdy.A21();
            H(npar_mapping + 1, 1) = -dypdy.A22();
            indices[npar_mapping] = fs->getIndexInMatrix();
            indices.at(npar_mapping + 1) = fs->getIndexInMatrix() + 1;
            ipar += npar_pos;
        }
        /* only consider proper motions of objects allowed to move,
        unless the fit is going to be degenerate */
        // TODO: left as reference for when we implement PM fitting
        // TODO: mjd here would become either deltaYears or maybe associations.epoch? Check the math first!
        // if (_fittingPM && fs->mightMove) {
        //     H(ipar, 0) = -mjd;  // Sign unchecked but consistent with above
        //     H(ipar + 1, 1) = -mjd;
        //     indices[ipar] = fs->getIndexInMatrix() + 2;
        //     indices[ipar + 1] = fs->getIndexInMatrix() + 3;
        //     ipar += npar_pm;
        // }

        // We can now compute the residual
        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);

        // do not write grad = H*transW*res to avoid
        // dynamic allocation of a temporary
        halpha = H * alpha;
        HW = H * transW;
        grad = HW * res;
        // now feed in triplets and fullGrad
        for (std::size_t ipar = 0; ipar < npar_tot; ++ipar) {
            for (std::size_t ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tripletList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            fullGrad(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }                    // end loop on measurements
    tripletList.setNextFreeIndex(kTriplets);
}

void AstrometryFit::leastSquareDerivativesReference(FittedStarList const &fittedStarList,
                                                    TripletList &tripletList,
                                                    Eigen::VectorXd &fullGrad) const {
    /**********************************************************************/
    /* @note the math in this method and accumulateStatRefStars() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* We compute here the derivatives of the terms involving fitted
       stars and reference stars. They only provide contributions if we
       are fitting positions: */
    if (!_fittingPos) return;
    /* the other case where the accumulation of derivatives stops
       here is when there are no RefStars */
    if (_associations->refStarList.size() == 0) return;
    Eigen::Matrix2d W(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::Matrix2d H(2, 2), halpha(2, 2), HW(2, 2);
    AstrometryTransformLinear der;
    Eigen::Vector2d res, grad;
    Eigen::Index indices[2];
    Eigen::Index kTriplets = tripletList.getNextFreeIndex();
    /* We cannot use the spherical coordinates directly to evaluate
       Euclidean distances, we have to use a projector on some plane in
       order to express least squares. Not projecting could lead to a
       disaster around the poles or across alpha=0.  So we need a
       projector. We construct a projector and will change its
       projection point at every object */
    TanRaDecToPixel proj(AstrometryTransformLinear(), Point(0., 0.));
    for (auto const &i : fittedStarList) {
        const FittedStar &fs = *i;
        auto rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        // Compute the derivative of the projector to incorporate its effects on the errors.
        proj.computeDerivative(fs, der, 1e-4);
        // sign checked. TODO check that the off-diagonal terms are OK.
        H(0, 0) = -der.A11();
        H(1, 0) = -der.A12();
        H(0, 1) = -der.A21();
        H(1, 1) = -der.A22();
        // TO DO : account for proper motions.
        double det = rsProj.vx * rsProj.vy - std::pow(rsProj.vxy, 2);
        if (rsProj.vx <= 0 || rsProj.vy <= 0 || det <= 0) {
            LOGLS_WARN(_log, "RefStar error matrix not positive definite for:  " << *rs);
            continue;
        }
        W(0, 0) = rsProj.vy / det;
        W(0, 1) = W(1, 0) = -rsProj.vxy / det;
        W(1, 1) = rsProj.vx / det;
        // compute alpha, a triangular square root
        // of W (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(W(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = W(0, 1) / alpha(0, 0);
        alpha(1, 1) = 1. / sqrt(det * W(0, 0));
        alpha(0, 1) = 0;
        indices[0] = fs.getIndexInMatrix();
        indices[1] = fs.getIndexInMatrix() + 1;
        unsigned npar_tot = 2;
        /* TODO: account here for proper motions in the reference
        catalog. We can code the effect and set the value to 0. Most
        (all?)  catalogs do not even come with a reference epoch. Gaia
        will change that. When refraction enters into the game, one should
        pay attention to the orientation of the frame */

        /* The residual should be Proj(fs)-Proj(*rs) in order to be consistent
        with the measurement terms. Since P(fs) = 0, we have: */
        res[0] = -rsProj.x;
        res[1] = -rsProj.y;
        halpha = H * alpha;
        // grad = H*W*res
        HW = H * W;
        grad = HW * res;
        // now feed in triplets and fullGrad
        for (std::size_t ipar = 0; ipar < npar_tot; ++ipar) {
            for (unsigned ic = 0; ic < 2; ++ic) {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
                tripletList.addTriplet(indices[ipar], kTriplets + ic, val);
            }
            fullGrad(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2;  // each measurement contributes 2 columns in the Jacobian
    }
    tripletList.setNextFreeIndex(kTriplets);
}

void AstrometryFit::accumulateStatImage(CcdImage const &ccdImage, Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesMeasurement() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/
    /* Setup */
    // 1 : get the Mapping's
    const AstrometryMapping *mapping = _astrometryModel->getMapping(ccdImage);
    // FittedStar is "observed" epoch, MeasuredStar is "baseline"
    double deltaYears = _epoch - ccdImage.getEpoch();
    // transformation from sky to TP
    auto sky2TP = _astrometryModel->getSkyToTangentPlane(ccdImage);
    // reserve matrix once for all measurements
    Eigen::Matrix2Xd transW(2, 2);

    auto &catalog = ccdImage.getCatalogForFit();
    for (auto const &ms : catalog) {
        if (!ms->isValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = *ms;
        tweakAstromMeasurementErrors(inPos, *ms, _posError);

        FatPoint outPos;
        // should *not* fill H if whatToFit excludes mapping parameters.
        mapping->transformPosAndErrors(inPos, outPos);
        double det = outPos.vx * outPos.vy - std::pow(outPos.vxy, 2);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, " Inconsistent measurement errors :drop measurement at "
                                     << Point(*ms) << " in image " << ccdImage.getName());
            continue;
        }
        transW(0, 0) = outPos.vy / det;
        transW(1, 1) = outPos.vx / det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy / det;

        std::shared_ptr<FittedStar const> const fs = ms->getFittedStar();
        Point fittedStarInTP = transformFittedStar(*fs, *sky2TP, deltaYears);

        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);
        double chi2Val = res.transpose() * transW * res;

        accum.addEntry(chi2Val, 2, ms);
    }  // end of loop on measurements
}

void AstrometryFit::accumulateStatImageList(CcdImageList const &ccdImageList, Chi2Accumulator &accum) const {
    for (auto const &ccdImage : ccdImageList) {
        accumulateStatImage(*ccdImage, accum);
    }
}

void AstrometryFit::accumulateStatRefStars(Chi2Accumulator &accum) const {
    /**********************************************************************/
    /** @note the math in this method and leastSquareDerivativesReference() must be kept consistent,
     * in terms of +/- convention, definition of model, etc. */
    /**********************************************************************/

    /* If you wonder why we project here, read comments in
       AstrometryFit::leastSquareDerivativesReference(TripletList &TList, Eigen::VectorXd &Rhs) */
    FittedStarList &fittedStarList = _associations->fittedStarList;
    TanRaDecToPixel proj(AstrometryTransformLinear(), Point(0., 0.));
    for (auto const &fs : fittedStarList) {
        auto rs = fs->getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(*fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        // TO DO : account for proper motions.
        double chi2 = computeProjectedRefStarChi2(rsProj);
        accum.addEntry(chi2, 2, fs);
    }
}

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void AstrometryFit::getIndicesOfMeasuredStar(MeasuredStar const &measuredStar, IndexVector &indices) const {
    if (_fittingDistortions) {
        const AstrometryMapping *mapping = _astrometryModel->getMapping(measuredStar.getCcdImage());
        mapping->getMappingIndices(indices);
    }
    std::shared_ptr<FittedStar const> const fs = measuredStar.getFittedStar();
    Eigen::Index fsIndex = fs->getIndexInMatrix();
    if (_fittingPos) {
        indices.push_back(fsIndex);
        indices.push_back(fsIndex + 1);
    }
    // For securing the outlier removal, the next block is just useless
    if (_fittingPM) {
        for (std::size_t k = 0; k < 2; ++k) indices.push_back(fsIndex + 2 + k);
    }
    /* Should not put the index of refaction stuff or we will not be
       able to remove more than 1 star at a time. */
}

void AstrometryFit::assignIndices(std::string const &whatToFit) {
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: Now fitting " << whatToFit);
    _fittingDistortions = (_whatToFit.find("Distortions") != std::string::npos);
    _fittingPos = (_whatToFit.find("Positions") != std::string::npos);
    _fittingPM = (_whatToFit.find("PM") != std::string::npos);
    // When entering here, we assume that whatToFit has already been interpreted.

    _nModelParams = 0;
    if (_fittingDistortions) _nModelParams = _astrometryModel->assignIndices(_whatToFit, 0);
    std::size_t ipar = _nModelParams;

    if (_fittingPos) {
        FittedStarList &fittedStarList = _associations->fittedStarList;
        for (auto &fittedStar : fittedStarList) {
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in GetMeasuredStarIndices
            fittedStar->setIndexInMatrix(ipar);
            ipar += 2;
            // TODO: left as reference for when we implement PM fitting
            // if ((_fittingPM)&fittedStar->mightMove) ipar += 2;
        }
    }
    _nStarParams = ipar - _nModelParams;
    _nTotal = ipar;
    LOGLS_DEBUG(_log, "nParameters total: " << _nTotal << " model: " << _nModelParams
                                            << " values: " << _nStarParams);
}

void AstrometryFit::offsetParams(Eigen::VectorXd const &delta) {
    if (delta.size() != _nTotal)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "AstrometryFit::offsetParams : the provided vector length is not compatible with "
                          "the current whatToFit setting");
    if (_fittingDistortions) _astrometryModel->offsetParams(delta);

    if (_fittingPos) {
        FittedStarList &fittedStarList = _associations->fittedStarList;
        for (auto const &i : fittedStarList) {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when assigning indices (assignIndices())
            Eigen::Index index = fs.getIndexInMatrix();
            fs.x += delta(index);
            fs.y += delta(index + 1);
            // TODO: left as reference for when we implement PM fitting
            // if ((_fittingPM)&fs.mightMove) {
            //     fs.pmx += delta(index + 2);
            //     fs.pmy += delta(index + 3);
            // }
        }
    }
}

void AstrometryFit::checkStuff() {
    const char *what2fit[] = {"Positions", "Distortions", "Positions Distortions"};
    // DEBUG
    for (unsigned k = 0; k < sizeof(what2fit) / sizeof(what2fit[0]); ++k) {
        assignIndices(what2fit[k]);
        TripletList tripletList(10000);
        Eigen::VectorXd grad(_nTotal);
        grad.setZero();
        leastSquareDerivatives(tripletList, grad);
        SparseMatrixD jacobian(_nTotal, tripletList.getNextFreeIndex());
        jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
        SparseMatrixD hessian = jacobian * jacobian.transpose();
        LOGLS_DEBUG(_log, "npar : " << _nTotal << ' ' << _nModelParams);
    }
}

void AstrometryFit::saveChi2MeasContributions(std::string const &filename) const {
    std::ofstream ofile(filename.c_str());
    std::string separator = "\t";

    ofile << "id" << separator << "xccd" << separator << "yccd " << separator;
    ofile << "rx" << separator << "ry" << separator;
    ofile << "xtp" << separator << "ytp" << separator;
    ofile << "mag" << separator << "deltaYears" << separator;
    ofile << "xErr" << separator << "yErr" << separator << "xyCov" << separator;
    ofile << "xtpi" << separator << "ytpi" << separator;
    ofile << "rxi" << separator << "ryi" << separator;
    ofile << "fsindex" << separator;
    ofile << "ra" << separator << "dec" << separator;
    ofile << "chi2" << separator << "nm" << separator;
    ofile << "detector" << separator << "visit" << std::endl;

    ofile << "id in source catalog" << separator << "coordinates in CCD (pixels)" << separator << separator;
    ofile << "residual on TP (degrees)" << separator << separator;
    ofile << "transformed coordinate in TP (degrees)" << separator << separator;
    ofile << "rough magnitude" << separator << "Julian epoch year delta from fit epoch" << separator;
    ofile << "transformed measurement uncertainty (degrees)" << separator << separator << separator;
    ofile << "as-read position on TP (degrees)" << separator << separator;
    ofile << "as-read residual on TP (degrees)" << separator << separator;
    ofile << "unique index of the fittedStar" << separator;
    ofile << "on sky position of fittedStar" << separator << separator;
    ofile << "contribution to chi2 (2D dofs)" << separator << "number of measurements of this fittedStar"
          << separator;
    ofile << "detector id" << separator << "visit id" << std::endl;

    const CcdImageList &ccdImageList = _associations->getCcdImageList();
    for (auto const &ccdImage : ccdImageList) {
        const MeasuredStarList &cat = ccdImage->getCatalogForFit();
        const AstrometryMapping *mapping = _astrometryModel->getMapping(*ccdImage);
        const auto readTransform = ccdImage->getReadWcs();
        // FittedStar is "observed" epoch, MeasuredStar is "baseline"
        double deltaYears = _epoch - ccdImage->getEpoch();
        for (auto const &ms : cat) {
            if (!ms->isValid()) continue;
            FatPoint tpPos;
            FatPoint inPos = *ms;
            tweakAstromMeasurementErrors(inPos, *ms, _posError);
            mapping->transformPosAndErrors(inPos, tpPos);
            auto sky2TP = _astrometryModel->getSkyToTangentPlane(*ccdImage);
            const std::unique_ptr<AstrometryTransform> readPixToTangentPlane =
                    compose(*sky2TP, *readTransform);
            FatPoint inputTpPos = readPixToTangentPlane->apply(inPos);
            std::shared_ptr<FittedStar const> const fs = ms->getFittedStar();

            Point fittedStarInTP = transformFittedStar(*fs, *sky2TP, deltaYears);
            Point res = tpPos - fittedStarInTP;
            Point inputRes = inputTpPos - fittedStarInTP;
            double det = tpPos.vx * tpPos.vy - std::pow(tpPos.vxy, 2);
            double wxx = tpPos.vy / det;
            double wyy = tpPos.vx / det;
            double wxy = -tpPos.vxy / det;
            double chi2 = wxx * res.x * res.x + wyy * res.y * res.y + 2 * wxy * res.x * res.y;

            ofile << std::setprecision(9);
            ofile << ms->getId() << separator << ms->x << separator << ms->y << separator;
            ofile << res.x << separator << res.y << separator;
            ofile << tpPos.x << separator << tpPos.y << separator;
            ofile << fs->getMag() << separator << deltaYears << separator;
            ofile << tpPos.vx << separator << tpPos.vy << separator << tpPos.vxy << separator;
            ofile << inputTpPos.x << separator << inputTpPos.y << separator;
            ofile << inputRes.x << separator << inputRes.y << separator;
            ofile << fs->getIndexInMatrix() << separator;
            ofile << fs->x << separator << fs->y << separator;
            ofile << chi2 << separator << fs->getMeasurementCount() << separator;
            ofile << ccdImage->getCcdId() << separator << ccdImage->getVisit() << std::endl;
        }  // loop on measurements in image
    }      // loop on images
}

void AstrometryFit::saveChi2RefContributions(std::string const &filename) const {
    std::ofstream ofile(filename.c_str());
    std::string separator = "\t";

    ofile << "ra" << separator << "dec " << separator;
    ofile << "rx" << separator << "ry" << separator;
    ofile << "mag" << separator;
    ofile << "xErr" << separator << "yErr" << separator << "xyCov" << separator;
    ofile << "fsindex" << separator;
    ofile << "chi2" << separator << "nm" << std::endl;

    ofile << "coordinates of fittedStar (degrees)" << separator << separator;
    ofile << "residual on TP (degrees)" << separator << separator;
    ofile << "magnitude" << separator;
    ofile << "refStar transformed measurement uncertainty (degrees)" << separator << separator << separator;
    ofile << "unique index of the fittedStar" << separator;
    ofile << "refStar contribution to chi2 (2D dofs)" << separator
          << "number of measurements of this FittedStar" << std::endl;

    // The following loop is heavily inspired from AstrometryFit::computeChi2()
    const FittedStarList &fittedStarList = _associations->fittedStarList;
    TanRaDecToPixel proj(AstrometryTransformLinear(), Point(0., 0.));
    for (auto const &i : fittedStarList) {
        const FittedStar &fs = *i;
        auto rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.setTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.transformPosAndErrors(*rs, rsProj);
        double chi2 = computeProjectedRefStarChi2(rsProj);

        ofile << std::setprecision(9);
        ofile << fs.x << separator << fs.y << separator;
        ofile << rsProj.x << separator << rsProj.y << separator;
        ofile << fs.getMag() << separator;
        ofile << rsProj.vx << separator << rsProj.vy << separator << rsProj.vxy << separator;
        ofile << fs.getIndexInMatrix() << separator;
        ofile << chi2 << separator << fs.getMeasurementCount() << std::endl;
    }  // loop on FittedStars
}
}  // namespace jointcal
}  // namespace lsst
