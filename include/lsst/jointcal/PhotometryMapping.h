// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
#define LSST_JOINTCAL_PHOTOMETRY_MAPPING_H

#include <memory>

#include "lsst/afw/image/PhotoCalib.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

/**
 * Relates transfo(s) to their position in the fitting matrix and allows interaction with the transfo(s).
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
    virtual unsigned getNpar() const = 0;

    /**
     * Return the on-sky transformed flux for measuredStar on ccdImage.
     *
     * @param[in]  measuredStar  The measured star position to transform.
     * @param[in]  instFlux      The instrument flux to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
     */
    virtual double transform(MeasuredStar const &measuredStar, double instFlux) const = 0;

    /**
     * Return the on-sky transformed flux uncertainty for measuredStar on ccdImage.
     * Matches the underlying PhotometryTransfo's `transformError()` until freezeErrorTransform() is called.
     *
     * @param[in]  measuredStar  The measured star position to transform.
     * @param[in]  value  The value to transform.
     * @param[in]  valueErr  The value uncertainty to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
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
     * @param[in]  instFlux     The instrument flux to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                             Eigen::Ref<Eigen::VectorXd> derivatives) const = 0;

    /**
     * Offset the transfo parameters by delta.
     *
     * @param[in]   delta vector to offset transfo parameters. Same ordering as derivatives in
     *              computeParameterDerivatives();
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    /// Make this mapping's parameters fixed (i.e. not varied during fitting).
    void setFixed(bool _fixed) { fixed = _fixed; }
    bool isFixed() { return fixed; }

    virtual Eigen::VectorXd getParameters() = 0;

    /**
     * Gets how this set of parameters (of length getNpar()) map into the "grand" fit.
     * Expects that indices has enough space reserved.
     */
    virtual void getMappingIndices(std::vector<unsigned> &indices) const = 0;

    /// Dump the contents of the transfos, for debugging.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /// Get the index of this mapping in the grand fit.
    unsigned getIndex() { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(unsigned i) { index = i; }

protected:
    // Start index of this mapping in the "grand" fit
    unsigned index;
    // Should this mapping be varied during fitting?
    bool fixed;
};

/**
 * A mapping containing a single photometryTransfo.
 */
class PhotometryMapping : public PhotometryMappingBase {
public:
    /**
     * Value transform takes ownership of transfo, error transform aliases it.
     *
     * Call freezeErrorTransform() to unalias the error transform.
     */
    explicit PhotometryMapping(std::shared_ptr<PhotometryTransfo> transfo)
            : PhotometryMappingBase(), _transfo(std::move(transfo)), _transfoErrors(_transfo) {}

    /// @copydoc PhotometryMappingBase::getNpar
    unsigned getNpar() const override {
        if (fixed) {
            return 0;
        } else {
            return _transfo->getNpar();
        }
    }

    /// @copydoc PhotometryMappingBase::transform
    double transform(MeasuredStar const &measuredStar, double value) const override {
        return _transfo->transform(measuredStar.x, measuredStar.y, value);
    }

    /// @copydoc PhotometryMappingBase::transformError
    double transformError(MeasuredStar const &measuredStar, double value, double valueErr) const override {
        return _transfoErrors->transformError(measuredStar.x, measuredStar.y, value, valueErr);
    }

    /// @copydoc PhotometryMappingBase::freezeErrorTransform
    void freezeErrorTransform() override {
        _transfoErrors = std::shared_ptr<PhotometryTransfo>(_transfo->clone());
    }

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double value,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        if (fixed) {
            return;
        } else {
            _transfo->computeParameterDerivatives(measuredStar.x, measuredStar.y, value, derivatives);
        }
    }

    /// @copydoc PhotometryMappingBase::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override { _transfo->offsetParams(delta); }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override { return _transfo->getParameters(); }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(std::vector<unsigned> &indices) const override {
        if (indices.size() < getNpar()) indices.resize(getNpar());
        for (unsigned k = 0; k < getNpar(); ++k) {
            indices[k] = index + k;
        }
    }

    /// @copydoc PhotometryMappingBase::dump
    void dump(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << " fixed: " << fixed << " transfo parameters: ";
        _transfo->dump(stream);
    }

    std::shared_ptr<PhotometryTransfo> getTransfo() const { return _transfo; }

    std::shared_ptr<PhotometryTransfo> getTransfoErrors() const { return _transfoErrors; }

private:
    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransfo> _transfo;
    // the transformation used for errors
    std::shared_ptr<PhotometryTransfo> _transfoErrors;
};

/**
 * A two-level photometric transform: one for the ccd and one for the visit.
 *
 * NOTE: We store the error in this mapping for now, because we are just using the single frame processing
 * calibration uncertainty, which is defined per Sensor Exposure, thus we would not preserve each
 * such uncertainty if we stored it in the chipMapping or visitMapping transfos individually.
 */
class ChipVisitPhotometryMapping : public PhotometryMappingBase {
public:
    ChipVisitPhotometryMapping(std::shared_ptr<PhotometryMapping> chipMapping,
                               std::shared_ptr<PhotometryMapping> visitMapping, double err)
            : PhotometryMappingBase(),
              _chipMapping(std::move(chipMapping)),
              _visitMapping(std::move(visitMapping)),
              _err(err) {}

    /// @copydoc PhotometryMappingBase::getNpar
    unsigned getNpar() const override { return _chipMapping->getNpar() + _visitMapping->getNpar(); }

    /// @copydoc PhotometryMappingBase::transform
    double transform(MeasuredStar const &measuredStar, double instFlux) const override {
        double tempFlux = _chipMapping->getTransfo()->transform(measuredStar.x, measuredStar.y, instFlux);
        return _visitMapping->getTransfo()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(),
                                                      tempFlux);
    }

    /// @copydoc PhotometryMappingBase::transformError
    double transformError(MeasuredStar const &measuredStar, double instFlux,
                          double instFluxErr) const override;

    /// @copydoc PhotometryMappingBase::freezeErrorTransform
    void freezeErrorTransform() override {
        _chipMapping->freezeErrorTransform();
        _visitMapping->freezeErrorTransform();
    }

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override;

    /// @copydoc PhotometryMappingBase::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override {
        _chipMapping->offsetParams(delta.segment(0, _chipMapping->getNpar()));
        _visitMapping->offsetParams(delta.segment(_chipMapping->getNpar(), _visitMapping->getNpar()));
    }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override {
        Eigen::VectorXd joined(getNpar());
        joined << _chipMapping->getParameters(), _visitMapping->getParameters();
        return joined;
    }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(std::vector<unsigned> &indices) const override;

    /// @copydoc PhotometryMappingBase::dump
    void dump(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << " chipMapping: ";
        _chipMapping->dump(stream);
        stream << "visitMapping: ";
        _visitMapping->dump(stream);
    }

    std::shared_ptr<PhotometryMapping> getChipMapping() const { return _chipMapping; }
    std::shared_ptr<PhotometryMapping> getVisitMapping() const { return _visitMapping; }

private:
    // the actual transformation to be fit
    std::shared_ptr<PhotometryMapping> _chipMapping;
    std::shared_ptr<PhotometryMapping> _visitMapping;

    // The uncertainty on this Exposure's calibration, taken from single frame processing.
    double _err;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
