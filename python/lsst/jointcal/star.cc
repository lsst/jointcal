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

#include <list>

#include "lsst/utils/python.h"

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

void declarePoint(py::module &mod) {
    py::class_<Point, std::shared_ptr<Point>> cls(mod, "Point");

    cls.def(py::init<double, double>(), "x"_a, "y"_a);

    cls.def_readwrite("x", &Point::x);
    cls.def_readwrite("y", &Point::y);

    utils::python::addOutputOp(cls, "__str__");
}

void declareFatPoint(py::module &mod) {
    py::class_<FatPoint, std::shared_ptr<FatPoint>, Point> cls(mod, "FatPoint");

    cls.def(py::init());
}

void declareBaseStar(py::module &mod) {
    py::class_<BaseStar, std::shared_ptr<BaseStar>, FatPoint> cls(mod, "BaseStar");

    cls.def(py::init<double, double, double, double>(), "x"_a, "y"_a, "flux"_a, "fluxErr"_a);

    // these three are actually declared in FatPoint, but we didn't want that in Python.
    // NOTE: see DM-9814 about the necessity of the pointer cast below.
    // readwrite so that we can set them in unittests.
    cls.def_readwrite("vx", (double BaseStar::*)&BaseStar::vx);
    cls.def_readwrite("vy", (double BaseStar::*)&BaseStar::vy);
    cls.def_readwrite("vxy", (double BaseStar::*)&BaseStar::vxy);

    // cls.def("getFlux", &BaseStar::getFlux);
    cls.def_property_readonly("flux", (double (BaseStar::*)() const) & BaseStar::getFlux);
    cls.def_property_readonly("fluxErr", (double (BaseStar::*)() const) & BaseStar::getFluxErr);
    cls.def_property_readonly("mag", (double (BaseStar::*)() const) & BaseStar::getMag);
    cls.def_property_readonly("magErr", (double (BaseStar::*)() const) & BaseStar::getMagErr);
}

void declareRefStar(py::module &mod) {
    py::class_<RefStar, std::shared_ptr<RefStar>, BaseStar> cls(mod, "RefStar");

    cls.def(py::init<double, double, double, double>(), "xx"_a, "yy"_a, "flux"_a, "fluxErr"_a);

    cls.def("setProperMotion", py::overload_cast<ProperMotion const &>(&RefStar::setProperMotion));
    cls.def("applyProperMotion", &RefStar::applyProperMotion);
}

void declareFittedStar(py::module &mod) {
    py::class_<FittedStar, std::shared_ptr<FittedStar>, BaseStar> cls(mod, "FittedStar");

    cls.def(py::init<>());
    cls.def(py::init<BaseStar const &>(), "baseStar"_a);
}

void declareMeasuredStar(py::module &mod) {
    py::class_<MeasuredStar, std::shared_ptr<MeasuredStar>, BaseStar> cls(mod, "MeasuredStar");

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
}

void declareProperMotion(py::module &mod) {
    py::class_<ProperMotion, std::shared_ptr<ProperMotion>> cls(mod, "ProperMotion");

    cls.def(py::init<double, double, double, double, double>(), "ra"_a, "dec"_a, "raErr"_a, "decErr"_a,
            "raDecCov"_a = 0);

    cls.def("apply", &ProperMotion::apply);

    utils::python::addOutputOp(cls, "__str__");
}

PYBIND11_MODULE(star, mod) {
    declarePoint(mod);
    declareFatPoint(mod);
    declareBaseStar(mod);
    declareRefStar(mod);
    declareFittedStar(mod);
    declareMeasuredStar(mod);
    declareProperMotion(mod);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
