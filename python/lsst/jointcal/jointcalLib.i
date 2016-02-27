// -*- lsst-c++ -*-
%define meas_jointcalLib_DOCSTRING
"
Python interface to lsst::jointcal classes
"
%enddef

%feature("autodoc", "1");

%module(package="lsst.jointcal", docstring=meas_jointcalLib_DOCSTRING) jointcalLib

%{
#include <exception>
#include <list>
#include <boost/shared_ptr.hpp>
#include "lsst/afw/table.h"
#include "lsst/afw/detection.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/jointcal/test.h"
#include "lsst/jointcal/test2.h"
#include "lsst/jointcal/jointcal.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstromFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/DistortionModel.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/jointcal/CountedRef.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/PhotomFit.h"
#include "lsst/jointcal/SimplePhotomModel.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_jointcal)

%import "lsst/afw/table/tableLib.i"

%include "lsst/pex/config.h"

%shared_ptr(lsst::daf::base::PropertySet);
%shared_ptr(lsst::jointcal::JointcalControl);

%template(SourceList) std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> >;
%template(WcsList) std::vector<boost::shared_ptr<lsst::afw::image::TanWcs> >;
%template(PropertySetList) std::vector<boost::shared_ptr<lsst::daf::base::PropertySet> >;
%template(CalibList) std::vector<boost::shared_ptr< lsst::afw::image::Calib > >;
%template(BboxList) std::vector<lsst::afw::geom::Box2I>;
%template(StringList) std::vector<std::string>;
%template(IntList) std::vector<int>;

%include "lsst/jointcal/test.h"
%include "lsst/jointcal/test2.h"
%include "lsst/jointcal/jointcal.h"
namespace lsst {
namespace jointcal {
class DistortionModel;
class TripletList;
class CcdImage;
class CcdImageList;
class MeasuredStar;
class FittedStarList;
class MeasuredStarList;
}}
%include "lsst/jointcal/Chi2.h"
%include "lsst/jointcal/AstromFit.h"
namespace lsst {
namespace jointcal {
class RefStarList;
class FittedStarList;
class Point;
}}
%include "lsst/jointcal/Associations.h"
namespace lsst {
namespace jointcal {
  class Gtransfo;
  class Mapping;
class SimplePolyMapping;
 class ProjectionHandler;
}}
%include "lsst/jointcal/DistortionModel.h"
namespace lsst {
namespace jointcal {
  class TanSipPix2RaDec;
}}
%include "lsst/jointcal/Projectionhandler.h"

namespace lsst {
namespace jointcal {
  class MeasuredStarList;
  class Frame;
}}

%shared_ptr(lsst::jointcal::CcdImage);
%template(CcdImageRef) boost::shared_ptr<lsst::jointcal::CcdImage>;
%template(CcdImageRefList) std::list<boost::shared_ptr<lsst::jointcal::CcdImage> >;
%template(TanSipPix2RaDecRef) boost::shared_ptr<lsst::jointcal::TanSipPix2RaDec>;

%shared_ptr(lsst::jointcal::TanSipPix2RaDec);
%include "lsst/jointcal/Point.h"
%include "lsst/jointcal/FatPoint.h"
%include "lsst/jointcal/Gtransfo.h"

%include "lsst/jointcal/CcdImage.h"
%include "lsst/jointcal/SimplePolyModel.h"
%include "lsst/jointcal/ConstrainedPolyModel.h"

%include "lsst/jointcal/SipToGtransfo.h"
%include "lsst/jointcal/PhotomModel.h"
%include "lsst/jointcal/PhotomFit.h"
%include "lsst/jointcal/SimplePhotomModel.h"

