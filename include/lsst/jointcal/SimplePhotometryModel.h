// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H

#include <map>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;

//! Photometric response model which has a single photometric factor per CcdImage.
class SimplePhotometryModel : public PhotometryModel {
public:
    SimplePhotometryModel(CcdImageList const &ccdImageList) { _myMap.reserve(ccdImageList.size()); }

    /// No copy or move: there is only ever one instance of a given model.
    SimplePhotometryModel(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel(SimplePhotometryModel &&) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel &&) = delete;

    /// @copydoc PhotometryModel::assignIndices
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) override;

    /// @copydoc PhotometryModel::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryModel::freezeErrorTransform
    void freezeErrorTransform() override;

    /// @copydoc PhotometryModel::getMappingIndices
    void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) const override;

    /// @copydoc PhotometryModel::getTotalParameters
    int getTotalParameters() const override;

    /// @copydoc PhotometryModel::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) const override;

    /// @copydoc PhotometryModel::dump
    void dump(std::ostream &stream = std::cout) const override;

protected:
    typedef std::unordered_map<CcdImageKey, std::unique_ptr<PhotometryMapping>> MapType;
    MapType _myMap;

    /// Return the mapping associated with this ccdImage.
    PhotometryMappingBase *findMapping(CcdImage const &ccdImage) const override;
};

class SimpleFluxModel : public SimplePhotometryModel {
public:
    SimpleFluxModel(CcdImageList const &ccdImageList);

    /// @copydoc PhotometryModel::offsetFittedStar
    void offsetFittedStar(FittedStar &fittedStar, double delta) const override {
        fittedStar.getFlux() -= delta;
    }

    /// @copydoc PhotometryModel::computeResidual
    double computeResidual(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transform
    double transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transformError
    double transformError(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::getRefError
    double getRefError(RefStar const &refStar) const override { return refStar.getFluxErr(); }

    /// @copydoc PhotometryModel::computeRefResidual
    double computeRefResidual(FittedStar const &fittedStar, RefStar const &refStar) const override {
        return fittedStar.getFlux() - refStar.getFlux();
    };

    /**
     * @copydoc PhotometryModel::toPhotoCalib
     *
     * @note SimplePhotometryModel uses a spatially-invariant transfo, so we can simplify the PhotoCalib.
     */
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const override;
};

class SimpleMagnitudeModel : public SimplePhotometryModel {
public:
    SimpleMagnitudeModel(CcdImageList const &ccdImageList);

    /// @copydoc PhotometryModel::offsetFittedStar
    void offsetFittedStar(FittedStar &fittedStar, double delta) const override {
        fittedStar.getMag() -= delta;
    }

    /// @copydoc PhotometryModel::computeResidual
    double computeResidual(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transform
    double transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transformError
    double transformError(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::getRefError
    double getRefError(RefStar const &refStar) const override { return refStar.getMagErr(); }

    /// @copydoc PhotometryModel::computeRefResidual
    double computeRefResidual(FittedStar const &fittedStar, RefStar const &refStar) const override {
        return fittedStar.getMag() - refStar.getMag();
    };

    /**
     * @copydoc PhotometryModel::toPhotoCalib
     *
     * @note SimplePhotometryModel uses a spatially-invariant transfo, so we can simplify the PhotoCalib.
     */
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const override;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
