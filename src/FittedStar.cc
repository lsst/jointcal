#include <iostream>
#include <iomanip>

#include "lsst/log/Log.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/StarList.h"

namespace {
LOG_LOGGER log = LOG_GET("jointcal.FastFinder");
}

namespace lsst {
namespace jointcal {

// cannot be in fittedstar.h, because of "crossed includes"
FittedStar::FittedStar(const MeasuredStar &measuredStar)
        : BaseStar(measuredStar),
          _mag(measuredStar.getMag()),
          _gen(-1),
          _wmag(measuredStar.getMagWeight()),
          _indexInMatrix(-1),
          _measurementCount(0),
          _refStar(nullptr) {
    _fluxErr = measuredStar.getInstFluxErr();
}

void FittedStar::setRefStar(const RefStar *refStar) {
    if ((_refStar != nullptr) && (refStar != nullptr)) {
        // TODO: should we raise an Exception in this case?
        LOGLS_WARN(log,
                   "FittedStar: " << *this << " is already matched to another RefStar. Clean up your lists.");
        LOGLS_WARN(log, "old refStar: " << *_refStar);
        LOGLS_WARN(log, "new refStar: " << *refStar);
    } else
        _refStar = refStar;
}

void FittedStar::addMagMeasurement(double magValue, double magWeight) {
    _mag = (_mag * _wmag + magValue * magWeight) / (_wmag + magWeight);
    _wmag += magWeight;
}

/************* FittedStarList ************************/

BaseStarList &fitted2Base(FittedStarList &starList) { return (BaseStarList &)starList; }

BaseStarList *fitted2Base(FittedStarList *starList) { return (BaseStarList *)starList; }

const BaseStarList &fitted2Base(const FittedStarList &starList) { return (const BaseStarList &)starList; }

const BaseStarList *fitted2Base(const FittedStarList *starList) { return (BaseStarList *)starList; }
}  // namespace jointcal
}  // namespace lsst
