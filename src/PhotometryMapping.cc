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
    // TODO DM-12161: possible optimization is to merge transform and computeDerivatives,
    // and/or save these intermediate calculations when transforming flux to use in derivatives.
    // Like what AstrometryMappings do with `computeTransformAndDerivatives` vs. `transformPosAndErrors`.

    double chipScale = _chipMapping->getTransfo()->transform(measuredStar.x, measuredStar.y, 1);
    double visitScale =
            _visitMapping->getTransfo()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(), 1);

    // NOTE: chipBlock is the product of the chip derivatives and the visit transforms, and vice versa.
    // NOTE: See DMTN-036 for the math behind this.
    if (_nParChips > 0 && !_chipMapping->isFixed()) {
        // The chip derivatives start at 0, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> chipBlock = derivatives.segment(0, _nParChips);
        _chipMapping->getTransfo()->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux,
                                                                chipBlock);
        chipBlock *= visitScale;
    }
    if (_nParVisits > 0) {
        // The visit derivatives start at the last chip derivative, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> visitBlock = derivatives.segment(_nParChips, _nParVisits);
        _visitMapping->getTransfo()->computeParameterDerivatives(
                measuredStar.getXFocal(), measuredStar.getYFocal(), instFlux, visitBlock);
        visitBlock *= chipScale;
    }
}

void ChipVisitPhotometryMapping::getMappingIndices(std::vector<unsigned> &indices) const {
    if (indices.size() < getNpar()) indices.resize(getNpar());
    if (_nParChips > 0) {
        _chipMapping->getMappingIndices(indices);
    }
    if (_nParVisits > 0) {
        // TODO DM-12169: there is probably a better way to feed a subpart of a std::vector
        // (maybe a view or iterators?)
        std::vector<unsigned> tempIndices(_visitMapping->getNpar());
        _visitMapping->getMappingIndices(tempIndices);
        for (unsigned k = 0; k < _visitMapping->getNpar(); ++k) {
            indices.at(k + _nParChips) = tempIndices.at(k);
        }
    }
}

void ChipVisitPhotometryMapping::setWhatToFit(bool const fittingChips, bool const fittingVisits) {
    if (fittingChips) {
        _nParChips = _chipMapping->getNpar();
    } else {
        _nParChips = 0;
    }
    if (fittingVisits) {
        _nParVisits = _visitMapping->getNpar();
    } else {
        _nParVisits = 0;
    }
}

}  // namespace jointcal
}  // namespace lsst
