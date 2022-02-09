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
    SimplePhotometryModel(CcdImageList const &ccdImageList, LOG_LOGGER log, double errorPedestal = 0)
            : PhotometryModel(log, errorPedestal) {
        _myMap.reserve(ccdImageList.size());
    }

    /// No copy or move: there is only ever one instance of a given model.
    SimplePhotometryModel(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel(SimplePhotometryModel &&) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel &&) = delete;

    /// @copydoc PhotometryModel::assignIndices
    Eigen::Index assignIndices(std::string const &whatToFit, Eigen::Index firstIndex) override;

    /// @copydoc PhotometryModel::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryModel::freezeErrorTransform
    void freezeErrorTransform() override;

    /// @copydoc PhotometryModel::getMappingIndices
    void getMappingIndices(CcdImage const &ccdImage, IndexVector &indices) const override;

    /// @copydoc PhotometryModel::getTotalParameters
    std::size_t getTotalParameters() const override;

    /// @copydoc PhotometryModel::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) const override;

    /// @copydoc PhotometryModel::print
    virtual void print(std::ostream &out) const override;

    ~SimplePhotometryModel() = default;

protected:
    using MapType = std::unordered_map<CcdImageKey, std::unique_ptr<PhotometryMapping>>;
    MapType _myMap;

    /// Return the mapping associated with this ccdImage.
    PhotometryMappingBase *findMapping(CcdImage const &ccdImage) const override;
};

class SimpleFluxModel : public SimplePhotometryModel {
public:
    SimpleFluxModel(CcdImageList const &ccdImageList, double errorPedestal = 0);

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
     * @note SimplePhotometryModel uses a spatially-invariant transform, so we can simplify the PhotoCalib.
     */
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const override;

    /// @copydoc PhotometryModel::print
    void print(std::ostream &out) const override;
};

class SimpleMagnitudeModel : public SimplePhotometryModel {
public:
    SimpleMagnitudeModel(CcdImageList const &ccdImageList, double errorPedestal = 0);

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
     * @note SimplePhotometryModel uses a spatially-invariant transform, so we can simplify the PhotoCalib.
     */
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const override;

    /// @copydoc PhotometryModel::print
    void print(std::ostream &out) const override;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
