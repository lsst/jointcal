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

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/AstrometryFit.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Chi2.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryFit(py::module &mod) {
    py::class_<AstrometryFit, std::shared_ptr<AstrometryFit>> cls(mod, "AstrometryFit");

    cls.def(py::init<Associations &, AstrometryModel *, double>(), "associations"_a, "astrometryModel"_a,
            "posError"_a);

    cls.def("minimize", &AstrometryFit::minimize, "whatToFit"_a, "nSigRejCut"_a = 0);
    cls.def("computeChi2", &AstrometryFit::computeChi2);
    cls.def("makeResTuple", &AstrometryFit::makeResTuple);
}

PYBIND11_PLUGIN(astrometryFit) {
    py::module::import("lsst.jointcal.associations");
    py::module::import("lsst.jointcal.astrometryModels");
    py::module::import("lsst.jointcal.chi2");
    py::module mod("astrometryFit");

    declareAstrometryFit(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
