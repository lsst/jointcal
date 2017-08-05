#include "lsst/log/Log.h"

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.PhotometryMapping");
}

namespace lsst {
namespace jointcal {

void ChipVisitPhotometryMapping::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                             double instFlux,
                                                             Eigen::Ref<Eigen::VectorXd> derivatives) const {
    // TODO: possible optimization is to merge transform and computeDerivatives,
    // and/or save these intermediate calculations when transforming flux to use in derivatives.
    // Like what AstrometryMappings do with `comptueTransformAndDerivatives` vs. `transformPosAndErrors`.

    double chipScale = _chipMapping->getTransfo()->transform(measuredStar.x, measuredStar.y, 1);
    double visitScale =
            _visitMapping->getTransfo()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(), 1);

    // NOTE: the chip derivatives start at 0, and the visit derivatives start at the last chip derivative.
    // This is independent of the full-fit indices.
    Eigen::Ref<Eigen::VectorXd> chipBlock = derivatives.segment(0, _chipMapping->getNpar());
    Eigen::Ref<Eigen::VectorXd> visitBlock =
            derivatives.segment(_chipMapping->getNpar(), _visitMapping->getNpar());

    if (not _chipMapping->isFixed()) {
        _chipMapping->getTransfo()->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux,
                                                                chipBlock);
        chipBlock *= visitScale;
    }
    _visitMapping->getTransfo()->computeParameterDerivatives(measuredStar.getXFocal(),
                                                             measuredStar.getYFocal(), instFlux, visitBlock);
    visitBlock *= chipScale;
}

void ChipVisitPhotometryMapping::getMappingIndices(std::vector<unsigned> &indices) const {
    if (indices.size() < getNpar()) indices.resize(getNpar());
    _chipMapping->getMappingIndices(indices);
    // TODO: there is probably a more elegant way to feed a subpart of a std::vector (can I get a view?)
    std::vector<unsigned> tempIndices(_visitMapping->getNpar());
    _visitMapping->getMappingIndices(tempIndices);
    for (unsigned k = 0; k < _visitMapping->getNpar(); ++k) {
        indices.at(k + _chipMapping->getNpar()) = tempIndices.at(k);
    }
}

}  // namespace jointcal
}  // namespace lsst
