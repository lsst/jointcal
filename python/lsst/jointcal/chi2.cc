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

#include "lsst/cpputils/python.h"

#include "lsst/jointcal/Chi2.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareChi2(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyChi2Statistic = py::class_<Chi2Statistic, std::shared_ptr<Chi2Statistic>>;

    wrappers.wrapType(PyChi2Statistic(wrappers.module, "Chi2Statistic"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cpputils::python::addOutputOp(cls, "__str__");
        cpputils::python::addOutputOp(cls, "__repr__");
        cls.def_readwrite("chi2", &Chi2Statistic::chi2);
        cls.def_readwrite("ndof", &Chi2Statistic::ndof);
    });
}
}  // namespace

void wrapChi2(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareChi2(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
