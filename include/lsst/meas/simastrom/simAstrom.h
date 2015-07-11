// -*- lsst-c++ -*-
#if !defined(LSST_MEAS_SIMASTROM_SIMASTROM_H)
#define LSST_MEAS_SIMASTROM_SIMASTROM_H

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace meas {
namespace simastrom {
    
class simAstrom {
public:
    
    simAstrom (
        std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > const sourceList,
        std::vector<PTR(lsst::daf::base::PropertySet)> const metaList,
        std::vector<PTR(lsst::afw::image::TanWcs)> const wcsList,
        std::vector<lsst::afw::geom::Box2I> const bboxList,
        std::vector<std::string> const filterList,
        std::vector<PTR(lsst::afw::image::Calib)> const calibList
    );
    
private:
    
    std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > _sourceList;
    std::vector <boost::shared_ptr<lsst::daf::base::PropertySet> > _metaList;
    std::vector<PTR(lsst::afw::image::TanWcs)> _wcsList;
    std::vector<lsst::afw::geom::Box2I> _bboxList;
    std::vector<std::string> _filterList;
    std::vector<PTR(lsst::afw::image::Calib)> _calibList;
};
    
}}}

#endif
