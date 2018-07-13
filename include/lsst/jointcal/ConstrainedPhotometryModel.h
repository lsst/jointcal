// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H

#include <map>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace lsst {
namespace jointcal {

/**
 * Photometry model with constraints, @f$M(x,y) = M_CCD(x,y)*M_visit(u,v)@f$
 *
 * This model consists of the following components:
 *   - A spatially invariant zero point per CCD, constrained across all visits, @f$M_CCD@f$.
 *   - A Chebyshev polynomial ( @f$a_ij*T_i(x)*T_j(y)@f$ ) per visit, constrained across all CCDs,
 *     @f$M_visit@f$.
 *
 * Because this model's parameters are degenerate under multiplication by a constant,
 * @f$M=(a*M_CCD)*(1/a*M_visit)@f$, we hold one CCD's zero point fixed to remove that degeneracy.
 */
class ConstrainedPhotometryModel : public PhotometryModel {
public:
    /**
     * Construct a constrained photometry model.
     *
     * @param      ccdImageList   The list of CCDImages to construct the model for.
     * @param      focalPlaneBBox The bounding box of the camera's focal plane, defining the domain of the
     *                            visit polynomial.
     * @param[in]  visitOrder    The order of the visit polynomial.
     */
    explicit ConstrainedPhotometryModel(CcdImageList const &ccdImageList,
                                        afw::geom::Box2D const &focalPlaneBBox, int visitOrder = 7);

    /// No copy or move: there is only ever one instance of a given model (i.e. per ccd+visit)
    ConstrainedPhotometryModel(ConstrainedPhotometryModel const &) = delete;
    ConstrainedPhotometryModel(ConstrainedPhotometryModel &&) = delete;
    ConstrainedPhotometryModel &operator=(ConstrainedPhotometryModel const &) = delete;
    ConstrainedPhotometryModel &operator=(ConstrainedPhotometryModel &&) = delete;

    /// @copydoc PhotometryModel::assignIndices
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) override;

    /// @copydoc PhotometryModel::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryModel::computeResidual
    double computeResidual(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transform
    double transform(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::transformError
    double transformError(CcdImage const &ccdImage, MeasuredStar const &measuredStar) const override;

    /// @copydoc PhotometryModel::freezeErrorTransform
    void freezeErrorTransform() override;

    /// @copydoc PhotometryModel::getMappingIndices
    void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) const override;

    /// @copydoc PhotometryModel::getTotalParameters
    int getTotalParameters() const override;

    /// @copydoc PhotometryModel::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) const override;

    /// @copydoc PhotometryModel::toPhotoCalib
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib(CcdImage const &ccdImage) const override;

    /// @copydoc PhotometryModel::dump
    void dump(std::ostream &stream = std::cout) const override;

private:
    PhotometryMappingBase *findMapping(CcdImage const &ccdImage) const override;

    // Which components of the model are we fitting currently?
    bool _fittingChips;
    bool _fittingVisits;

    /* The per-ccdImage transforms, each of which is a composition of a chip and visit transform.
     * Not all pairs of _visitMap[visit] and _chipMap[chip] are guaranteed to have an entry in
     * _chipVisitMap (for example, one ccd in one visit might fail to produce a catalog).
     */
    typedef std::unordered_map<CcdImageKey, std::unique_ptr<ChipVisitPhotometryMapping>> MapType;
    MapType _chipVisitMap;

    // The per-visit transforms that go into _chipVisitMap.
    typedef std::map<VisitIdType, std::shared_ptr<PhotometryMapping>> VisitMapType;
    VisitMapType _visitMap;
    // The per-chip transforms that go into _chipVisitMap.
    typedef std::map<CcdIdType, std::shared_ptr<PhotometryMapping>> ChipMapType;
    ChipMapType _chipMap;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CONSTRAINED_PHOTOMETRY_MODEL_H
