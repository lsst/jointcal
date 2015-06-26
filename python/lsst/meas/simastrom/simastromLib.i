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
#include "lsst/daf/base/PropertySet.h"
#include "lsst/meas/simastrom/test.h"
#include "lsst/meas/simastrom/test2.h"
#include "lsst/meas/simastrom/simAstrom.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_simastrom)

%import "lsst/afw/table/tableLib.i"

%shared_ptr(lsst::daf::base::PropertySet);

%template(SourceList) std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> >;
%template(WcsList) std::vector<boost::shared_ptr<lsst::afw::image::Wcs> >;
%template(PropertySetList) std::vector<boost::shared_ptr<lsst::daf::base::PropertySet> >;

%include "lsst/meas/simastrom/test.h"
%include "lsst/meas/simastrom/test2.h"
%include "lsst/meas/simastrom/simAstrom.h"

