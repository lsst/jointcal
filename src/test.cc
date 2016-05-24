#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "boost/shared_array.hpp"
#include "boost/multi_index_container.hpp"
#include "boost/multi_index/sequenced_index.hpp"
#include "boost/multi_index/ordered_index.hpp"
#include "boost/multi_index/global_fun.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/jointcal/test.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace jointcal {
    
    void test(
        lsst::afw::table::SourceCatalog const &sourceCat,
        lsst::daf::base::PropertySet const &metaData
    ) {
        std::cout << "Source catalog size : " << sourceCat.size() << std::endl;
        double lat = metaData.get<double>("LATITUDE");
        std::cout << "Observatory Latitude : " << lat << std::endl;
        int count = 0;
        for (afwTable::SourceCatalog::const_iterator sourcePtr = sourceCat.begin();
            sourcePtr != sourceCat.end(); ++sourcePtr) {
                if (count < 2) {
                    double x = sourcePtr->getX();
                    std::cout << "Gaussian centroid X : " << x  << std::endl;
                }
                count++;
        }
    }
}} // end of namespaces
