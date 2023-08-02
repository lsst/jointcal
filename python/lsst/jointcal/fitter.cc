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

void declareFitterBase(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyFitterBase = py::class_<FitterBase>;

    wrappers.wrapType(PyFitterBase(wrappers.module, "FitterBase"), [](auto &mod, auto &cls) {

        cls.def("minimize", &FitterBase::minimize, "whatToFit"_a, "nSigRejCut"_a = 0,
                "sigmaRelativeTolerance"_a = 0, "doRankUpdate"_a = true, "doLineSearch"_a = false,
                "dumpMatrixFile"_a = "");
        cls.def("computeChi2", &FitterBase::computeChi2);
        cls.def("saveChi2Contributions", &FitterBase::saveChi2Contributions);
    });
}

void declareAstrometryFit(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryFit = py::class_<AstrometryFit, FitterBase>;

    wrappers.wrapType(PyAstrometryFit(wrappers.module, "AstrometryFit"), [](auto &mod, auto &cls) {
        cls.def(py::init<std::shared_ptr<Associations>, std::shared_ptr<AstrometryModel>, double>(),
                "associations"_a, "astrometryModel"_a, "posError"_a);
        cls.def("getModel", &AstrometryFit::getModel, py::return_value_policy::reference_internal);
    });
}

void declarePhotometryFit(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryFit = py::class_<PhotometryFit, FitterBase>;

    wrappers.wrapType(PyPhotometryFit(wrappers.module, "PhotometryFit"), [](auto &mod, auto &cls) {
        cls.def(py::init<std::shared_ptr<Associations>, std::shared_ptr<PhotometryModel>>(), "associations"_a,
                "photometryModel"_a);

        cls.def("getModel", &PhotometryFit::getModel, py::return_value_policy::reference_internal);
    });
}
}  // namespace

void wrapFitter(lsst::cpputils::python::WrapperCollection &wrappers) {
    wrappers.wrapType(py::enum_<MinimizeResult>(wrappers.module, "MinimizeResult"), [](auto &mod, auto &enm) {
        enm.value("Converged", MinimizeResult::Converged);
        enm.value("Chi2Increased", MinimizeResult::Chi2Increased);
        enm.value("NonFinite", MinimizeResult::NonFinite);
        enm.value("Failed", MinimizeResult::Failed);
    });

    declareFitterBase(wrappers);
    declareAstrometryFit(wrappers);
    declarePhotometryFit(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
