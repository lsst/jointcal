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

    // NOTE: the chip derivatives start at 0, and the visit derivatives start at the last chip derivative.
    // This is independent of the full-fit indices.
    Eigen::Ref<Eigen::VectorXd> chipBlock = derivatives.segment(0, _chipMapping->getNpar());
    Eigen::Ref<Eigen::VectorXd> visitBlock =
            derivatives.segment(_chipMapping->getNpar(), _visitMapping->getNpar());

    // NOTE: chipBlock is the product of the chip derivatives and the visit transforms, and vice versa.
    // NOTE: See DMTN-036 for the math behind this.
    if (!_chipMapping->isFixed()) {
        _chipMapping->getTransfo()->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux,
                                                                chipBlock);
        chipBlock *= visitScale;
    }
    _visitMapping->getTransfo()->computeParameterDerivatives(measuredStar.getXFocal(),
                                                             measuredStar.getYFocal(), instFlux, visitBlock);
    visitBlock *= chipScale;
}

double ChipVisitPhotometryMapping::transformError(MeasuredStar const &measuredStar, double instFlux,
                                                  double instFluxErr) const {
    // NOTE: we assume here that the only source of error from the mapping itself is the uncertainty from
    // single frame processing (i.e. the uncertainty on f_0^-1).

    double tempFlux = _chipMapping->getTransfoErrors()->transform(measuredStar.x, measuredStar.y, instFlux);
    double flux = _visitMapping->getTransfoErrors()->transform(measuredStar.getXFocal(),
                                                               measuredStar.getYFocal(), tempFlux);
    auto chipErrTransfo = _chipMapping->getTransfoErrors();
    auto visitErrTransfo = _visitMapping->getTransfoErrors();
    // The calibration scale factor at the location of this source is the "true" scale factor that goes
    // into the uncertainty calculation. See PhotoCalib.h, which is defined as: flux = instFlux*scale(x,y).
    tempFlux = chipErrTransfo->transform(measuredStar.x, measuredStar.y, 1.0);
    double scale = visitErrTransfo->transform(measuredStar.getXFocal(), measuredStar.getYFocal(), tempFlux);

    return flux * std::hypot(_err / scale, instFluxErr / instFlux);
}

void ChipVisitPhotometryMapping::getMappingIndices(std::vector<unsigned> &indices) const {
    if (indices.size() < getNpar()) indices.resize(getNpar());
    _chipMapping->getMappingIndices(indices);
    // TODO DM-12169: there is probably a better way to feed a subpart of a std::vector (a view or iterators?)
    std::vector<unsigned> tempIndices(_visitMapping->getNpar());
    _visitMapping->getMappingIndices(tempIndices);
    for (unsigned k = 0; k < _visitMapping->getNpar(); ++k) {
        indices.at(k + _chipMapping->getNpar()) = tempIndices.at(k);
    }
}

}  // namespace jointcal
}  // namespace lsst
