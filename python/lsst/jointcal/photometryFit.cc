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
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/PhotometryFit.h"
#include "lsst/jointcal/PhotometryModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryFit(py::module &mod) {
    py::class_<PhotometryFit, std::shared_ptr<PhotometryFit>> cls(mod, "PhotometryFit");

    cls.def(py::init<Associations &, PhotometryModel *>(), "associations"_a, "photometryModel"_a);

    cls.def("minimize", &PhotometryFit::minimize, "whatToFit"_a);
    cls.def("computeChi2", &PhotometryFit::computeChi2);
    cls.def("makeResTuple", &PhotometryFit::makeResTuple);
}

PYBIND11_PLUGIN(photometryFit) {
    py::module::import("lsst.jointcal.associations");
    py::module::import("lsst.jointcal.chi2");
    py::module::import("lsst.jointcal.photometryModels");
    py::module mod("photometryFit");

    declarePhotometryFit(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
