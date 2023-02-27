// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "lsst/cpputils/python.h"

#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/RefStar.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePoint(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPoint = py::class_<Point, std::shared_ptr<Point>>;

    wrappers.wrapType(PyPoint(wrappers.module, "Point"), [](auto &mod, auto &cls) {
        cls.def(py::init<double, double>(), "x"_a, "y"_a);
        cls.def_readwrite("x", &Point::x);
        cls.def_readwrite("y", &Point::y);
        utils::python::addOutputOp(cls, "__str__");
    });
}

void declareFatPoint(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyFatPoint = py::class_<FatPoint, std::shared_ptr<FatPoint>, Point> ;

    wrappers.wrapType(PyFatPoint(wrappers.module, "FatPoint"), [](auto &mod, auto &cls) {
        cls.def(py::init());
    });
}

void declareBaseStar(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyBaseStar = py::class_<BaseStar, std::shared_ptr<BaseStar>, FatPoint>;

    wrappers.wrapType(PyBaseStar(wrappers.module, "BaseStar"), [](auto &mod, auto &cls) {
        cls.def(py::init<double, double, double, double>(), "x"_a, "y"_a, "flux"_a, "fluxErr"_a);

        // these three are actually declared in FatPoint, but we didn't want that in Python.
        // readwrite so that we can set them in unittests.
        cls.def_readwrite("vx", &BaseStar::vx);
        cls.def_readwrite("vy", &BaseStar::vy);
        cls.def_readwrite("vxy", &BaseStar::vxy);

        cls.def_property_readonly("flux", (double (BaseStar::*)() const) &BaseStar::getFlux);
        cls.def_property_readonly("fluxErr", (double (BaseStar::*)() const) &BaseStar::getFluxErr);
        cls.def_property_readonly("mag", (double (BaseStar::*)() const) &BaseStar::getMag);
        cls.def_property_readonly("magErr", (double (BaseStar::*)() const) &BaseStar::getMagErr);
    });
}

void declareRefStar(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyRefStar = py::class_<RefStar, std::shared_ptr<RefStar>, BaseStar>;

    wrappers.wrapType(PyRefStar(wrappers.module, "RefStar"), [](auto &mod, auto &cls) {
        cls.def(py::init<double, double, double, double>(), "xx"_a, "yy"_a, "flux"_a, "fluxErr"_a);
        cls.def("setProperMotion", py::overload_cast<ProperMotion const &>(&RefStar::setProperMotion));
        cls.def("applyProperMotion", &RefStar::applyProperMotion);
    });
}

void declareFittedStar(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyFittedStar = py::class_<FittedStar, std::shared_ptr<FittedStar>, BaseStar> ;

    wrappers.wrapType(PyFittedStar(wrappers.module, "FittedStar"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<BaseStar const &>(), "baseStar"_a);
    });
}

void declareMeasuredStar(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyMeasuredStar = py::class_<MeasuredStar, std::shared_ptr<MeasuredStar>, BaseStar>;

    wrappers.wrapType(PyMeasuredStar(wrappers.module, "MeasuredStar"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<BaseStar const &>(), "baseStar"_a);
        cls.def("getFittedStar", &MeasuredStar::getFittedStar);
        cls.def("setFittedStar", &MeasuredStar::setFittedStar);
        cls.def("getInstFlux", &MeasuredStar::getInstFlux);
        cls.def("setInstFluxAndErr", &MeasuredStar::setInstFluxAndErr);
        cls.def("getInstFluxErr", &MeasuredStar::getInstFluxErr);
        cls.def("getInstMag", &MeasuredStar::getInstMag);
        cls.def("getInstMagErr", &MeasuredStar::getInstMagErr);
        cls.def("setXFocal", &MeasuredStar::setXFocal);
        cls.def("setYFocal", &MeasuredStar::setYFocal);
        cls.def("getXFocal", &MeasuredStar::getXFocal);
        cls.def("getYFocal", &MeasuredStar::getYFocal);
    });
}

void declareProperMotion(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyProperMotion = py::class_<ProperMotion, std::shared_ptr<ProperMotion>> ;

    wrappers.wrapType(PyProperMotion(wrappers.module, "ProperMotion"), [](auto &mod, auto &cls) {
        cls.def(py::init<double, double, double, double, double>(), "ra"_a, "dec"_a, "raErr"_a, "decErr"_a,
                "raDecCov"_a = 0);
        cls.def("apply", &ProperMotion::apply);
        utils::python::addOutputOp(cls, "__str__");
    });
}
}  // namespace

void wrapStar(lsst::cpputils::python::WrapperCollection &wrappers) {
    declarePoint(wrappers);
    declareFatPoint(wrappers);
    declareBaseStar(wrappers);
    declareRefStar(wrappers);
    declareFittedStar(wrappers);
    declareMeasuredStar(wrappers);
    declareProperMotion(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
