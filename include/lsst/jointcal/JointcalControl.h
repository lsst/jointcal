// -*- lsst-c++ -*-
#if !defined(LSST_JOINTCAL_JointcalControl_H)
#define LSST_JOINTCAL_JointcalControl_H

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace jointcal {

struct JointcalControl {
    LSST_CONTROL_FIELD(sourceFluxField, std::string, "name of flux field in source catalog");

    explicit JointcalControl(std::string const& sourceFluxField = "slot_CalibFlux") :
        // Set sourceFluxType to the value used in the source selector.
        sourceFluxField(sourceFluxField)
    {
        validate();
    }

    ~JointcalControl() {};

    void validate() const {
        if (sourceFluxField.empty()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "sourceFluxField must be specified");
        }
    }
};

}
} // end of namespaces

#endif
