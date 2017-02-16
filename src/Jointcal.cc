#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "boost/make_shared.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/jointcal/Jointcal.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/AstromFit.h"

#include "Eigen/Core"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace jointcal {

void JointcalControl::validate() const {
    if (sourceFluxField.empty()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterError, "sourceFluxField must be specified");
    }
}

}
} // end of namespaces
