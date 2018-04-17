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
#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/SimpleAstrometryModel.h"
#include "lsst/jointcal/ConstrainedAstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryModel(py::module &mod) {
    py::class_<AstrometryModel, std::shared_ptr<AstrometryModel>> cls(mod, "AstrometryModel");

    cls.def("getNpar", &AstrometryModel::getNpar);
    cls.def("getMapping", &AstrometryModel::getMapping, py::return_value_policy::reference_internal);
    cls.def("assignIndices", &AstrometryModel::assignIndices);
    cls.def("offsetParams", &AstrometryModel::offsetParams);
    cls.def("getSky2TP", &AstrometryModel::getSky2TP);
    cls.def("makeSkyWcs", &AstrometryModel::makeSkyWcs);
}

void declareSimpleAstrometryModel(py::module &mod) {
    py::class_<SimpleAstrometryModel, std::shared_ptr<SimpleAstrometryModel>, AstrometryModel> cls(
            mod, "SimpleAstrometryModel");

    cls.def(py::init<CcdImageList const &, const std::shared_ptr<ProjectionHandler const>, bool, unsigned,
                     unsigned>(),
            "ccdImageList"_a, "projectionHandler"_a, "initFromWcs"_a, "nNotFit"_a = 0, "order"_a = 3);

    cls.def("getTransfo", &SimpleAstrometryModel::getTransfo, py::return_value_policy::reference_internal);
}

void declareConstrainedAstrometryModel(py::module &mod) {
    py::class_<ConstrainedAstrometryModel, std::shared_ptr<ConstrainedAstrometryModel>, AstrometryModel> cls(
            mod, "ConstrainedAstrometryModel");

    cls.def(py::init<CcdImageList const &, std::shared_ptr<ProjectionHandler const>, int, int>(),
            "ccdImageList"_a, "projectionHandler"_a, "chipOrder"_a, "visitOrder"_a);

    cls.def("getChipTransfo", &ConstrainedAstrometryModel::getChipTransfo,
            py::return_value_policy::reference_internal);
    cls.def("getVisitTransfo", &ConstrainedAstrometryModel::getVisitTransfo,
            py::return_value_policy::reference_internal);
}

PYBIND11_PLUGIN(astrometryModels) {
    py::module::import("lsst.jointcal.ccdImage");
    py::module::import("lsst.jointcal.gtransfo");
    py::module::import("lsst.jointcal.astrometryMappings");
    py::module mod("astrometryModels");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declareAstrometryModel(mod);
    declareSimpleAstrometryModel(mod);
    declareConstrainedAstrometryModel(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
