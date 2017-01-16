// -*- lsst-c++ -*-
#if !defined(LSST_JOINTCAL_JOINTCAL_H)
#define LSST_JOINTCAL_JOINTCAL_H

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/daf/base/PropertySet.h"

namespace lsst {
namespace jointcal {

struct JointcalControl {
    LSST_CONTROL_FIELD(sourceFluxField, std::string, "name of flux field in source catalog");

    JointcalControl(std::string const & sourceFluxField = "slot_CalibFlux") :
        // Set sourceFluxType to the value used in the source selector.
        sourceFluxField(sourceFluxField)
    {
        validate();
    }

    void validate() const;

    ~JointcalControl() {};
};

class Jointcal {
public:

    Jointcal (
        std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > const sourceList,
//        std::vector<PTR(lsst::afw::image::VisitInfo)> const visitInfoList,
        std::vector<PTR(lsst::daf::base::PropertySet)> const metaList,
        std::vector<PTR(lsst::afw::image::TanWcs)> const wcsList,
        std::vector<lsst::afw::geom::Box2I> const bboxList,
        std::vector<std::string> const filterList,
        std::vector<PTR(lsst::afw::image::Calib)> const calibList,
        std::vector<int> const visitList,
        std::vector<int> const ccdList,
        PTR(lsst::jointcal::JointcalControl) const control
    );

private:

    std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > _sourceList;
//    std::vector <std::shared_ptr<lsst::afw::image::VisitInfo> > _visitInfoList;
    std::vector <std::shared_ptr<lsst::daf::base::PropertySet> > _metaList;
    std::vector<PTR(lsst::afw::image::TanWcs)> _wcsList;
    std::vector<lsst::afw::geom::Box2I> _bboxList;
    std::vector<std::string> _filterList;
    std::vector<PTR(lsst::afw::image::Calib)> _calibList;
    std::vector<int> const _visitList;
    std::vector<int> const _ccdList;
};

}
} // end of namespaces

#endif
