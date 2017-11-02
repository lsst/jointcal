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
#include "lsst/jointcal/FitterBase.h"
#include "lsst/jointcal/PhotometryFit.h"
#include "lsst/jointcal/PhotometryModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareFitterBase(py::module &mod) {
    py::class_<FitterBase, std::shared_ptr<FitterBase>> cls(mod, "FitterBase");

    cls.def("minimize", &FitterBase::minimize, "whatToFit"_a, "nSigRejCut"_a = 0);
    cls.def("computeChi2", &FitterBase::computeChi2);
    cls.def("saveChi2Contributions", &FitterBase::saveChi2Contributions);
}

void declareAstrometryFit(py::module &mod) {
    py::class_<AstrometryFit, std::shared_ptr<AstrometryFit>, FitterBase> cls(mod, "AstrometryFit");

    cls.def(py::init<std::shared_ptr<Associations>, std::shared_ptr<AstrometryModel>, double>(),
            "associations"_a, "astrometryModel"_a, "posError"_a);
}

void declarePhotometryFit(py::module &mod) {
    py::class_<PhotometryFit, std::shared_ptr<PhotometryFit>, FitterBase> cls(mod, "PhotometryFit");

    cls.def(py::init<std::shared_ptr<Associations>, std::shared_ptr<PhotometryModel>>(), "associations"_a,
            "photometryModel"_a);
}

PYBIND11_PLUGIN(fitter) {
    py::module::import("lsst.jointcal.associations");
    py::module::import("lsst.jointcal.astrometryModels");
    py::module::import("lsst.jointcal.chi2");
    py::module::import("lsst.jointcal.photometryModels");
    py::module mod("fitter");

    py::enum_<MinimizeResult>(mod, "MinimizeResult")
            .value("Converged", MinimizeResult::Converged)
            .value("Chi2Increased", MinimizeResult::Chi2Increased)
            .value("Failed", MinimizeResult::Failed);

    declareFitterBase(mod);
    declareAstrometryFit(mod);
    declarePhotometryFit(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
