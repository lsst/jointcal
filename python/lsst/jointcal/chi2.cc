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

#include "lsst/jointcal/Chi2.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareChi2(py::module &mod) {
    py::class_<Chi2, std::shared_ptr<Chi2>> cls(mod, "Chi2");

    cls.def(py::init<>());

    cls.def("__str__", &Chi2::__str__);
    cls.def_readwrite("chi2", &Chi2::chi2);
    cls.def_readwrite("ndof", &Chi2::ndof);
}

PYBIND11_PLUGIN(chi2) {
    py::module mod("chi2");

    declareChi2(mod);

    return mod.ptr();
}

}}}  // lsst::jointcal::<anonymous>
