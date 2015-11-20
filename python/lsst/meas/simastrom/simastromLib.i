// -*- lsst-c++ -*-
%define meas_simastromLib_DOCSTRING
"
Python interface to lsst::meas::simastrom classes
"
%enddef

%feature("autodoc", "1");

%module(package="lsst.meas.simastrom", docstring=meas_simastromLib_DOCSTRING) simastromLib

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
#include "lsst/meas/simastrom/test.h"
#include "lsst/meas/simastrom/test2.h"
#include "lsst/meas/simastrom/simAstrom.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/AstromFit.h"
#include "lsst/meas/simastrom/Associations.h"
#include "lsst/meas/simastrom/DistortionModel.h"
#include "lsst/meas/simastrom/SimplePolyModel.h"
#include "lsst/meas/simastrom/Projectionhandler.h"
#include "lsst/meas/simastrom/CountedRef.h"
#include "lsst/meas/simastrom/SipToGtransfo.h"
#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/meas/simastrom/FatPoint.h"
#include "lsst/meas/simastrom/Point.h"
#include "lsst/meas/simastrom/ConstrainedPolyModel.h"
#include "lsst/meas/simastrom/PhotomFit.h"
#include "lsst/meas/simastrom/SimplePhotomModel.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_simastrom)

%import "lsst/afw/table/tableLib.i"

%include "lsst/pex/config.h"

%shared_ptr(lsst::daf::base::PropertySet);
%shared_ptr(lsst::meas::simastrom::SimAstromControl);

%template(SourceList) std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> >;
%template(WcsList) std::vector<boost::shared_ptr<lsst::afw::image::TanWcs> >;
%template(PropertySetList) std::vector<boost::shared_ptr<lsst::daf::base::PropertySet> >;
%template(CalibList) std::vector<boost::shared_ptr< lsst::afw::image::Calib > >;
%template(BboxList) std::vector<lsst::afw::geom::Box2I>;
%template(StringList) std::vector<std::string>;
%template(IntList) std::vector<int>;

%include "lsst/meas/simastrom/test.h"
%include "lsst/meas/simastrom/test2.h"
%include "lsst/meas/simastrom/simAstrom.h"
namespace lsst {
namespace meas {
namespace simastrom {
class DistortionModel;
class TripletList;
class CcdImage;
class CcdImageList;
class MeasuredStar;
class FittedStarList;
class MeasuredStarList;
}}}
%include "lsst/meas/simastrom/Chi2.h"
%include "lsst/meas/simastrom/AstromFit.h"
namespace lsst {
namespace meas {
namespace simastrom {
class RefStarList;
class FittedStarList;
class Point;
}}}
%include "lsst/meas/simastrom/Associations.h"
namespace lsst {
namespace meas {
namespace simastrom {
  class Gtransfo;
  class Mapping;
class SimplePolyMapping;
 class ProjectionHandler;
}}}
%include "lsst/meas/simastrom/DistortionModel.h"
namespace lsst {
namespace meas {
namespace simastrom {
  class TanSipPix2RaDec;
}}}
%include "lsst/meas/simastrom/Projectionhandler.h"

namespace lsst {
namespace meas {
namespace simastrom {
  class MeasuredStarList;
  class Frame;
}}}

%shared_ptr(lsst::meas::simastrom::CcdImage);
%template(CcdImageRef) boost::shared_ptr<lsst::meas::simastrom::CcdImage>;
%template(CcdImageRefList) std::list<boost::shared_ptr<lsst::meas::simastrom::CcdImage> >;
%template(TanSipPix2RaDecRef) boost::shared_ptr<lsst::meas::simastrom::TanSipPix2RaDec>;

%shared_ptr(lsst::meas::simastrom::TanSipPix2RaDec);
%include "lsst/meas/simastrom/Point.h"
%include "lsst/meas/simastrom/FatPoint.h"
%include "lsst/meas/simastrom/Gtransfo.h"

%include "lsst/meas/simastrom/CcdImage.h"
%include "lsst/meas/simastrom/SimplePolyModel.h"
%include "lsst/meas/simastrom/ConstrainedPolyModel.h"

%include "lsst/meas/simastrom/SipToGtransfo.h"
%include "lsst/meas/simastrom/PhotomModel.h"
%include "lsst/meas/simastrom/PhotomFit.h"
%include "lsst/meas/simastrom/SimplePhotomModel.h"

