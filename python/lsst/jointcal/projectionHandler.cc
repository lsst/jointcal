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

#include "lsst/utils/python.h"

#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/CcdImage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareProjectionHandler(py::module &mod) {
    py::class_<ProjectionHandler, std::shared_ptr<ProjectionHandler>> cls(mod, "ProjectionHandler");
    utils::python::addOutputOp(cls, "__str__");
    utils::python::addOutputOp(cls, "__repr__");
}

void declareIdentityProjectionHandler(py::module &mod) {
    py::class_<IdentityProjectionHandler, std::shared_ptr<IdentityProjectionHandler>, ProjectionHandler> cls(
            mod, "IdentityProjectionHandler");
    cls.def(py::init());
}

void declareOneTPPerVisitHandler(py::module &mod) {
    py::class_<OneTPPerVisitHandler, std::shared_ptr<OneTPPerVisitHandler>, ProjectionHandler> cls(
            mod, "OneTPPerVisitHandler");
    cls.def(py::init<CcdImageList const &>(), "ccdImageList"_a);
}

PYBIND11_MODULE(projectionHandler, mod) {
    py::module::import("lsst.jointcal.ccdImage");
    declareProjectionHandler(mod);
    declareIdentityProjectionHandler(mod);
    declareOneTPPerVisitHandler(mod);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
