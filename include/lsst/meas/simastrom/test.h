// -*- lsst-c++ -*-
#if !defined(LSST_MEAS_SIMASTROM_TEST_H)
#define LSST_MEAS_SIMASTROM_TEST_H

#include <cmath>
#include <string>
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"

namespace lsst {
namespace meas {
namespace simastrom {
    
    void test(
        lsst::afw::table::SourceCatalog const &sourceCat
    );
    
    
}}}

#endif
