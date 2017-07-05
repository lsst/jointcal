/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <list>

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

    cls.def_readonly("x", &Point::x);
    cls.def_readonly("y", &Point::y);
}

void declareBaseStar(py::module &mod) {
    py::class_<BaseStar, std::shared_ptr<BaseStar>, Point> cls(mod, "BaseStar");

    // these three are actually declared in FatPoint, but we don't need that in Python.
    // NOTE: see DM-9814 about the necessity of the pointer cast below.
    cls.def_readonly("vx", (double BaseStar::*)&BaseStar::vx);
    cls.def_readonly("vy", (double BaseStar::*)&BaseStar::vy);
    cls.def_readonly("vxy", (double BaseStar::*)&BaseStar::vxy);

    // cls.def("getFlux", &BaseStar::getFlux);
    cls.def_property_readonly("flux", (double (BaseStar::*)() const) & BaseStar::getFlux);

    cls.def("__str__", &BaseStar::__str__);
}

void declareRefStar(py::module &mod) {
    py::class_<RefStar, std::shared_ptr<RefStar>, BaseStar> cls(mod, "RefStar");

    cls.def(py::init<const BaseStar &>(), "baseStar"_a);
}

void declareFittedStar(py::module &mod) {
    py::class_<FittedStar, std::shared_ptr<FittedStar>, BaseStar> cls(mod, "FittedStar");

    cls.def(py::init<const BaseStar &>(), "baseStar"_a);
}

void declareMeasuredStar(py::module &mod) {
    py::class_<MeasuredStar, std::shared_ptr<MeasuredStar>, BaseStar> cls(mod, "MeasuredStar");

    cls.def(py::init<const BaseStar &>(), "baseStar"_a);
}

PYBIND11_PLUGIN(star) {
    py::module mod("star");

    declarePoint(mod);
    declareBaseStar(mod);
    declareRefStar(mod);
    declareFittedStar(mod);
    declareMeasuredStar(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
