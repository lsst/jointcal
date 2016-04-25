#ifndef LSST_MEAS_SIMASTROM_MATCHEXPOSURE_H
#define LSST_MEAS_SIMASTROM_MATCHEXPOSURE_H


#include "lsst/jointcal/ExposureCatalog.h"
#include "lsst/jointcal/ChipArrangement.h"

namespace lsst {
namespace jointcal {
    
  class ExposureCatalog;
  class Point;
  struct SimAstromControl;

  //! Routine to astrometrically match a whole exposure at once, relying on a ChipArrangement
  bool MatchExposure(ExposureCatalog &EC, const Point &TangentPoint, const JointcalControl &Control);
    
}}

#endif

