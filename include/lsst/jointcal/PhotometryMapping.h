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

#ifndef LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
#define LSST_JOINTCAL_PHOTOMETRY_MAPPING_H

#include <memory>

#include "lsst/afw/image/PhotoCalib.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransform.h"

namespace lsst {
namespace jointcal {

/**
 * Relates transform(s) to their position in the fitting matrix and allows interaction with the transform(s).
 */
class PhotometryMappingBase {
public:
    PhotometryMappingBase() : index(-1), fixed(false) {}
    virtual ~PhotometryMappingBase(){};

    /// No copy or move: there is only ever one instance of a given mapping (i.e. per ccd+visit)
    PhotometryMappingBase(PhotometryMappingBase const &) = delete;
    PhotometryMappingBase(PhotometryMappingBase &&) = delete;
    PhotometryMappingBase &operator=(PhotometryMappingBase const &) = delete;
    PhotometryMappingBase &operator=(PhotometryMappingBase &&) = delete;

    /// Number of total parameters in this mapping
    virtual std::size_t getNpar() const = 0;

    /**
     * Return the on-sky transformed flux for measuredStar on ccdImage.
     *
     * @param[in]  measuredStar  The measured star position to transform.
     * @param[in]  value         The instrument flux or magnitude to transform.
     *
     * @return     The on-sky value () transformed from value at measuredStar's position.
     */
    virtual double transform(MeasuredStar const &measuredStar, double value) const = 0;

    /**
     * Return the on-sky transformed flux uncertainty for measuredStar on ccdImage.
     * Matches the underlying PhotometryTransform's `transformError()` until freezeErrorTransform() is called.
     *
     * @param[in]  measuredStar  The measured star position to transform.
     * @param[in]  value  The flux or magnitude to transform.
     * @param[in]  valueErr  The flux or magnitude uncertainty to transform.
     *
     * @return     The on-sky value transformed from value at measuredStar's position.
     */
    virtual double transformError(MeasuredStar const &measuredStar, double value, double valueErr) const = 0;

    /**
     * Once this routine has been called, the error transform is not modified by offsetParams().
     *
     * The routine can be called when the mappings are roughly in place. After the call, the transformations
     * used to propagate errors are no longer affected when updating the mappings. This allows an exactly
     * linear fit, which can be necessary for some model+data combinations.
     */
    virtual void freezeErrorTransform() = 0;

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  measuredStar The measured star position to transform.
     * @param[in]  value        The instrument flux or magnitude to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void computeParameterDerivatives(MeasuredStar const &measuredStar, double value,
                                             Eigen::Ref<Eigen::VectorXd> derivatives) const = 0;

    /// Make this mapping's parameters fixed (i.e. not varied during fitting).
    void setFixed(bool _fixed) { fixed = _fixed; }
    bool isFixed() { return fixed; }

    virtual Eigen::VectorXd getParameters() = 0;

    /**
     * Gets how this set of parameters (of length getNpar()) map into the "grand" fit.
     * Expects that indices has enough space reserved.
     */
    virtual void getMappingIndices(IndexVector &indices) const = 0;

    /// Print the contents of the transforms, for debugging.
    virtual void print(std::ostream &stream = std::cout) const = 0;

    /// Get the index of this mapping in the grand fit.
    Eigen::Index getIndex() { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(Eigen::Index i) { index = i; }

protected:
    // Start index of this mapping in the "grand" fit
    Eigen::Index index;
    // Should this mapping be varied during fitting?
    bool fixed;
};

/**
 * A mapping containing a single photometryTransform.
 */
class PhotometryMapping : public PhotometryMappingBase {
public:
    /**
     * Value transform takes ownership of transform, error transform aliases it.
     *
     * Call freezeErrorTransform() to unalias the error transform.
     */
    explicit PhotometryMapping(std::shared_ptr<PhotometryTransform> transform)
            : PhotometryMappingBase(), _transform(std::move(transform)), _transformErrors(_transform) {}

    /// @copydoc PhotometryMappingBase::getNpar
    std::size_t getNpar() const override {
        if (fixed) {
            return 0;
        } else {
            return _transform->getNpar();
        }
    }

    /// @copydoc PhotometryMappingBase::transform
    double transform(MeasuredStar const &measuredStar, double value) const override {
        return _transform->transform(measuredStar.x, measuredStar.y, value);
    }

    /// @copydoc PhotometryMappingBase::transformError
    double transformError(MeasuredStar const &measuredStar, double value, double valueErr) const override {
        return _transformErrors->transformError(measuredStar.x, measuredStar.y, value, valueErr);
    }

    /// @copydoc PhotometryMappingBase::freezeErrorTransform
    void freezeErrorTransform() override {
        _transformErrors = std::shared_ptr<PhotometryTransform>(_transform->clone());
    }

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        if (fixed) {
            return;
        } else {
            _transform->computeParameterDerivatives(measuredStar.x, measuredStar.y, value, derivatives);
        }
    }

    /**
     * Offset the transform parameters by delta.
     *
     * @param[in]   delta vector to offset transform parameters. Same ordering as derivatives in
     *              computeParameterDerivatives();
     */
    void offsetParams(Eigen::VectorXd const &delta) { _transform->offsetParams(delta); }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override { return _transform->getParameters(); }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(IndexVector &indices) const override {
        if (indices.size() < getNpar()) indices.resize(getNpar());
        for (std::size_t k = 0; k < getNpar(); ++k) {
            indices[k] = index + k;
        }
    }

    /// @copydoc PhotometryMappingBase::print
    void print(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << " fixed: " << fixed << " transform parameters: ";
        _transform->print(stream);
    }

    std::shared_ptr<PhotometryTransform> getTransform() const { return _transform; }

    std::shared_ptr<PhotometryTransform> getTransformErrors() const { return _transformErrors; }

private:
    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransform> _transform;
    // the transformation used for errors
    std::shared_ptr<PhotometryTransform> _transformErrors;
};

/**
 * A two-level photometric transform: one for the ccd and one for the visit.
 */
class ChipVisitPhotometryMapping : public PhotometryMappingBase {
public:
    ChipVisitPhotometryMapping(std::shared_ptr<PhotometryMapping> chipMapping,
                               std::shared_ptr<PhotometryMapping> visitMapping)
            : PhotometryMappingBase(),
              _nParChip(0),
              _nParVisit(0),
              _chipMapping(std::move(chipMapping)),
              _visitMapping(std::move(visitMapping)) {}

    /// @copydoc PhotometryMappingBase::getNpar
    std::size_t getNpar() const override { return _nParChip + _nParVisit; }

    /// @copydoc PhotometryMappingBase::transform
    double transform(MeasuredStar const &measuredStar, double value) const override {
        double temp = _chipMapping->getTransform()->transform(measuredStar.x, measuredStar.y, value);
        return _visitMapping->getTransform()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(),
                                                        temp);
    }

    /// @copydoc PhotometryMappingBase::freezeErrorTransform
    void freezeErrorTransform() override {
        _chipMapping->freezeErrorTransform();
        _visitMapping->freezeErrorTransform();
    }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override {
        Eigen::VectorXd joined(getNpar());
        joined << _chipMapping->getParameters(), _visitMapping->getParameters();
        return joined;
    }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(IndexVector &indices) const override;

    /**
     * Set whether to fit chips or visits.
     *
     * This must be called before anything that depends on knowing the number of parameters in the fit,
     * such as offsetParams(), getParameters(), or computeParameterDerivatives().
     *
     * @param fittingChips Fit the chip transform.
     * @param fittingVisits Fit the visit transform.
     */
    void setWhatToFit(bool const fittingChips, bool const fittingVisits);

    /// @copydoc PhotometryMappingBase::print
    void print(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << std::endl << "chip mapping: ";
        _chipMapping->print(stream);
        stream << std::endl << "visit mapping: ";
        _visitMapping->print(stream);
    }

    std::shared_ptr<PhotometryMapping> getChipMapping() const { return _chipMapping; }
    std::shared_ptr<PhotometryMapping> getVisitMapping() const { return _visitMapping; }

    std::size_t getNParChip() const { return _nParChip; }
    std::size_t getNParVisit() const { return _nParVisit; }

protected:
    // These are either transform.getNpar() or 0, depending on whether we are fitting that component or not.
    std::size_t _nParChip, _nParVisit;

    // the actual transformation to be fit
    std::shared_ptr<PhotometryMapping> _chipMapping;
    std::shared_ptr<PhotometryMapping> _visitMapping;
};

class ChipVisitFluxMapping : public ChipVisitPhotometryMapping {
public:
    ChipVisitFluxMapping(std::shared_ptr<PhotometryMapping> chipMapping,
                         std::shared_ptr<PhotometryMapping> visitMapping)
            : ChipVisitPhotometryMapping(chipMapping, visitMapping) {}

    /// @copydoc PhotometryMappingBase::transformError
    double transformError(MeasuredStar const &measuredStar, double value, double valueErr) const override;

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override;
};

class ChipVisitMagnitudeMapping : public ChipVisitPhotometryMapping {
public:
    ChipVisitMagnitudeMapping(std::shared_ptr<PhotometryMapping> chipMapping,
                              std::shared_ptr<PhotometryMapping> visitMapping)
            : ChipVisitPhotometryMapping(chipMapping, visitMapping) {}

    /**
     * @copydoc PhotometryMappingBase::transformError
     *
     * @note This method takes instFlux and instFluxErr: the error calculation has to
     * use fluxes to get the math right.
     */
    double transformError(MeasuredStar const &measuredStar, double value, double valueErr) const override;

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
