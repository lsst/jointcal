#include <cmath>

#include "lsst/log/Log.h"

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.PhotometryMapping");
}

namespace lsst {
namespace jointcal {

void ChipVisitPhotometryMapping::getMappingIndices(std::vector<unsigned> &indices) const {
    if (indices.size() < getNpar()) indices.resize(getNpar());
    // If we're fitting the chip mapping, fill those indices.
    if (_nParChip > 0) {
        _chipMapping->getMappingIndices(indices);
    }
    // If we're fitting the visit mapping, fill those indices.
    if (_nParVisit > 0) {
        // TODO DM-12169: there is probably a better way to feed a subpart of a std::vector
        // (maybe a view or iterators?)
        std::vector<unsigned> tempIndices(_visitMapping->getNpar());
        _visitMapping->getMappingIndices(tempIndices);
        // We have to insert the visit indices starting after the chip indices.
        for (unsigned k = 0; k < _visitMapping->getNpar(); ++k) {
            indices.at(k + _nParChip) = tempIndices.at(k);
        }
    }
}

void ChipVisitPhotometryMapping::setWhatToFit(bool const fittingChips, bool const fittingVisits) {
    if (fittingChips) {
        _nParChip = _chipMapping->getNpar();
    } else {
        _nParChip = 0;
    }
    if (fittingVisits) {
        _nParVisit = _visitMapping->getNpar();
    } else {
        _nParVisit = 0;
    }
}

// ChipVisitFluxMapping methods

double ChipVisitFluxMapping::transformError(MeasuredStar const &measuredStar, double instFlux,
                                            double instFluxErr) const {
    // The transformed error is s_m = dM(f,x,y)/df + s_f.
    double tempFlux =
            _chipMapping->getTransfoErrors()->transform(measuredStar.x, measuredStar.y, instFluxErr);
    return _visitMapping->getTransfoErrors()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(),
                                                        tempFlux);
}

void ChipVisitFluxMapping::computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                                       Eigen::Ref<Eigen::VectorXd> derivatives) const {
    // TODO DM-12161: possible optimization is to merge transform and computeDerivatives,
    // and/or save these intermediate calculations when transforming flux to use in derivatives.
    // Like what AstrometryMappings do with `computeTransformAndDerivatives` vs. `transformPosAndErrors`.

    double chipScale = _chipMapping->getTransfo()->transform(measuredStar.x, measuredStar.y, 1);
    double visitScale =
            _visitMapping->getTransfo()->transform(measuredStar.getXFocal(), measuredStar.getYFocal(), 1);

    // NOTE: chipBlock is the product of the chip derivatives and the visit transforms, and vice versa.
    // NOTE: See DMTN-036 for the math behind this.
    if (getNParChip() > 0 && !_chipMapping->isFixed()) {
        // The chip derivatives start at 0, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> chipBlock = derivatives.segment(0, getNParChip());
        _chipMapping->getTransfo()->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux,
                                                                chipBlock);
        chipBlock *= visitScale;
    }
    if (getNParVisit() > 0) {
        // The visit derivatives start at the last chip derivative, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> visitBlock = derivatives.segment(getNParChip(), getNParVisit());
        _visitMapping->getTransfo()->computeParameterDerivatives(
                measuredStar.getXFocal(), measuredStar.getYFocal(), instFlux, visitBlock);
        visitBlock *= chipScale;
    }
}

// ChipVisitMagnitudeMapping methods

double ChipVisitMagnitudeMapping::transformError(MeasuredStar const &measuredStar, double instFlux,
                                                 double instFluxErr) const {
    // The transformed error is s_mout = 2.5/ln(10) * instFluxErr / instFlux
    // because the other components of the mapping (f0, the polynomials) disappear in the partial derivative.
    return 2.5 / std::log(10.0) * instFluxErr / instFlux;
}

void ChipVisitMagnitudeMapping::computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                                            Eigen::Ref<Eigen::VectorXd> derivatives) const {
    // TODO DM-12161: possible optimization is to merge transform and computeDerivatives,
    // and/or save these intermediate calculations when transforming flux to use in derivatives.
    // Like what AstrometryMappings do with `computeTransformAndDerivatives` vs. `transformPosAndErrors`.

    // NOTE: See DMTN-036 for the math behind this.
    if (getNParChip() > 0 && !_chipMapping->isFixed()) {
        // The chip derivatives start at 0, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> chipBlock = derivatives.segment(0, getNParChip());
        _chipMapping->getTransfo()->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux,
                                                                chipBlock);
    }
    if (getNParVisit() > 0) {
        // The visit derivatives start at the last chip derivative, independent of the full-fit indices.
        Eigen::Ref<Eigen::VectorXd> visitBlock = derivatives.segment(getNParChip(), getNParVisit());
        _visitMapping->getTransfo()->computeParameterDerivatives(
                measuredStar.getXFocal(), measuredStar.getYFocal(), instFlux, visitBlock);
    }
}

}  // namespace jointcal
}  // namespace lsst
