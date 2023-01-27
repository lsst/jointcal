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
 * This program is free software: you can redistribute it and/or wrappersify
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

#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/CcdImage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareProjectionHandler(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyProjectionHandler =  py::class_<ProjectionHandler, std::shared_ptr<ProjectionHandler>>;

    wrappers.wrapType(PyProjectionHandler(wrappers.module, "ProjectionHandler"), [](auto &mod, auto &cls) {
        utils::python::addOutputOp(cls, "__str__");
        utils::python::addOutputOp(cls, "__repr__");
    });
}

void declareIdentityProjectionHandler(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyIdentityProjectionHandler =
            py::class_<IdentityProjectionHandler, std::shared_ptr<IdentityProjectionHandler>, ProjectionHandler>;

    wrappers.wrapType(
            PyIdentityProjectionHandler(wrappers.module, "IdentityProjectionHandler"), [](auto &mod, auto &cls) {
                cls.def(py::init());
            });
}

void declareOneTPPerVisitHandler(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyOneTPPerVisitHandler=  py::class_<OneTPPerVisitHandler, std::shared_ptr<OneTPPerVisitHandler>, ProjectionHandler>;

    wrappers.wrapType(PyOneTPPerVisitHandler(wrappers.module, "OneTPPerVisitHandler"), [](auto &mod, auto &cls) {
        cls.def(py::init<CcdImageList const &>(), "ccdImageList"_a);
    });
}
}  // namespace

void wrapProjectionHandler(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareProjectionHandler(wrappers);
    declareIdentityProjectionHandler(wrappers);
    declareOneTPPerVisitHandler(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
