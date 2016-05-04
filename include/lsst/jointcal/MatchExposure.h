#ifndef LSST_MEAS_SIMASTROM_MATCHEXPOSURE_H
#define LSST_MEAS_SIMASTROM_MATCHEXPOSURE_H


#include "lsst/jointcal/ExposureCatalog.h"
#include "lsst/jointcal/ChipArrangement.h"
#include "lsst/jointcal/CountedRef.h"

namespace lsst {
namespace jointcal {
    
  class ExposureCatalog;
  class Point;
  struct SimAstromControl;

  // prefer this to a typedef because it saves a %template declaration to swig.
  struct ExposureSolutionType : public std::map<unsigned, lsst::jointcal::CountedRef<Gtransfo> >
    {
    };

  //! Routine to astrometrically match a whole exposure at once, relying on a ChipArrangement. The ourine returns mappings from pixel space to tangent plane.
ExposureSolutionType MatchExposure(ExposureCatalog &EC, const Point &TangentPoint, const JointcalControl &Control);
    
}}

#endif

