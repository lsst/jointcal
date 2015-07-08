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
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/meas/simastrom/simAstrom.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Point.h"

#include "Eigen/Core"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace meas {
namespace simastrom {
    
    simAstrom::simAstrom(
        std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > const sourceList,
        std::vector<PTR(lsst::daf::base::PropertySet)> const metaList,
        std::vector<PTR(lsst::afw::image::Wcs)> const wcsList,
        std::vector<lsst::afw::geom::Box2I> const bboxList
    ):
        _sourceList(sourceList),
        _metaList(metaList),
        _wcsList(wcsList),
        _bboxList(bboxList)
{
    std::cout << "simAstrom constructor invoked " << std::endl;
    std::cout << "Vectors contain : " << _sourceList.size() << " elements" << std::endl;
    std::cout << _metaList[0]->get<double>("LATITUDE") << std::endl;
    std::cout << _sourceList[1][10].getRa() << std::endl;
    std::cout << _wcsList[1]->getPixelOrigin() << std::endl;
    std::cout << _bboxList[1] << std::endl;
    
    // Check how to get SIP coefficients from WCS
    lsst::daf::base::PropertyList::Ptr wcsMeta = _wcsList[1]->getFitsMetadata();
//    std::cout << wcsMeta->getOrderedNames() << std::endl;
    std::cout << wcsMeta->get<int>("A_ORDER") << std::endl;
    
    Eigen::MatrixXd sipA;
    lsst:afw::image::TanWcs::decodeSipHeader(*wcsMeta, "A", sipA);
    std::cout << sipA << std::endl;
    
    Point where(100.,200.);
    CcdImage(_sourceList[1], where, _wcsList[1], _metaList[1], _bboxList[1]);
}    
    
}}}
