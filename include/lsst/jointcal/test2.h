// -*- lsst-c++ -*-
#if !defined(LSST_JOINTCAL_TEST2_H)
#define LSST_JOINTCAL_TEST2_H

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace jointcal {
    
class test2 {
public:
    
    test2 (
    );
    
    void addSourceTable (
        lsst::afw::table::SourceCatalog const &sourceCat,
        lsst::daf::base::PropertySet const &metaData
//        lsst::afw::image::Wcs const &wcs
    );

    void addMetaList(std::vector<PTR(lsst::daf::base::PropertySet)>);
    
private:
    
    typedef std::tuple<lsst::afw::table::SourceCatalog,
            lsst::daf::base::PropertySet> sourceTuple;
//    std::vector <lsst::afw::table::SourceCatalog> _sources;
//    std::vector <lsst::daf::base::PropertySet> _sources;
    std::vector <boost::shared_ptr<lsst::daf::base::PropertySet> > _sources;
};
    
}} // end of namespaces

#endif
