#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "boost/scoped_array.hpp"
#include "boost/shared_array.hpp"
#include "boost/multi_index_container.hpp"
#include "boost/multi_index/sequenced_index.hpp"
#include "boost/multi_index/ordered_index.hpp"
#include "boost/multi_index/global_fun.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/simastrom/test.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace meas {
namespace simastrom {
    
    void test(
        lsst::afw::table::SourceCatalog const &sourceCat
    ) {
        std::cout << "Source catalog size : " << sourceCat.size() << std::endl;
    }
}}}
