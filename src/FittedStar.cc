#include <iostream>
#include <iomanip>

#include "lsst/log/Log.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/StarList.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.FastFinder");
}

namespace lsst {
namespace jointcal {

// cannot be in fittedstar.h, because of "crossed includes"
FittedStar::FittedStar(const MeasuredStar &M)
        : BaseStar(M),
          mag(M.Mag()),
          emag(-1),
          col(0.),
          gen(-1),
          wmag(M.MagWeight()),
          indexInMatrix(-1),
          measurementCount(0),
          _refStar(nullptr),
          flux2(-1),
          fluxErr2(-1) {
    fluxErr = M.eflux;
}

void FittedStar::setRefStar(const RefStar *refStar) {
    if ((_refStar != nullptr) && (refStar != nullptr)) {
        // TODO: should we raise an Exception in this case?
        LOGLS_WARN(_log,
                   "FittedStar: " << *this << " is already matched to another RefStar. Clean up your lists.");
        LOGLS_WARN(_log, "old refStar: " << *_refStar);
        LOGLS_WARN(_log, "new refStar: " << *refStar);
    } else
        _refStar = refStar;
}

void FittedStar::AddMagMeasurement(double MagValue, double MagWeight) {
    mag = (mag * wmag + MagValue * MagWeight) / (wmag + MagWeight);
    wmag += MagWeight;
}

/************* FittedStarList ************************/

BaseStarList &Fitted2Base(FittedStarList &This) { return (BaseStarList &)This; }

BaseStarList *Fitted2Base(FittedStarList *This) { return (BaseStarList *)This; }

const BaseStarList &Fitted2Base(const FittedStarList &This) { return (const BaseStarList &)This; }

const BaseStarList *Fitted2Base(const FittedStarList *This) { return (BaseStarList *)This; }
}  // namespace jointcal
}  // namespace lsst
