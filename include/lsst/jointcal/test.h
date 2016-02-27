// -*- lsst-c++ -*-
#if !defined(LSST_JOINTCAL_TEST_H)
#define LSST_JOINTCAL_TEST_H

#include <cmath>
#include <string>
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace jointcal {
    
    void test(
        lsst::afw::table::SourceCatalog const &sourceCat,
        lsst::daf::base::PropertySet const &metaData
    );
    
    
}} // end of namespaces

#endif
