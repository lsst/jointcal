#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

//#include "boost/scoped_array.hpp"
//#include "boost/shared_array.hpp"
//#include "boost/multi_index_container.hpp"
//#include "boost/multi_index/sequenced_index.hpp"
//#include "boost/multi_index/ordered_index.hpp"
//#include "boost/multi_index/global_fun.hpp"

#include "boost/make_shared.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/jointcal/test2.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace jointcal {
    
    test2::test2()
{
    std::cout << "test2 constructor called " << std::endl;
}    
    
void test2::addSourceTable(
    lsst::afw::table::SourceCatalog const &sourceCat,
    lsst::daf::base::PropertySet const &metaData
//    lsst::afw::image::Wcs const &wcs
    ) {
//        _sources.push_back(sourceCat);
        PTR(lsst::daf::base::PropertySet) tt = metaData.deepCopy();
//        _sources.push_back(PTR(lsst::daf::base::PropertySet) metaData);
        _sources.push_back(tt);
//        _sources.push_back(metaData);
        std::cout << "tuple feeded with a new source table " << std::endl;
    }
    
void test2::addMetaList(std::vector<PTR(lsst::daf::base::PropertySet)> ll) {
    std::cout << ll[0]->get<double>("LATITUDE") << std::endl;
}
    
}} // end of namespaces
