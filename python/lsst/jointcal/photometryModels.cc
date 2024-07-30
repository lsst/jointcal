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
#include "pybind11/eigen.h"
#include "pybind11/stl.h"
#include "Eigen/Core"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/ConstrainedPhotometryModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryModel =  py::class_<PhotometryModel, std::shared_ptr<PhotometryModel>>;

    wrappers.wrapType(PyPhotometryModel(wrappers.module, "PhotometryModel"), [](auto &mod, auto &cls) {
        cls.def("assignIndices", &PhotometryModel::assignIndices);
        cls.def("freezeErrorTransform", &PhotometryModel::freezeErrorTransform);

        cls.def("offsetParams", &PhotometryModel::offsetParams);
        cls.def("offsetFittedStar", &PhotometryModel::offsetFittedStar);

        cls.def("transform", &PhotometryModel::transform);
        cls.def("transformError", &PhotometryModel::transformError);
        cls.def("computeResidual", &PhotometryModel::computeResidual);

        cls.def("getRefError", &PhotometryModel::getRefError);
        cls.def("computeRefResidual", &PhotometryModel::computeRefResidual);

        cls.def("checkPositiveOnBBox", &PhotometryModel::checkPositiveOnBBox);
        cls.def("validate", &PhotometryModel::validate);

        cls.def("getMappingIndices", &PhotometryModel::getMappingIndices);
        cls.def("computeParameterDerivatives",
                [](PhotometryModel const &self, MeasuredStar const &star, CcdImage const &ccdImage) {
                    Eigen::VectorXd derivatives(self.getNpar(ccdImage));
                    self.computeParameterDerivatives(star, ccdImage, derivatives);
                    return derivatives;
                });

        cls.def("getNpar", &PhotometryModel::getNpar);
        cls.def("toPhotoCalib", &PhotometryModel::toPhotoCalib);
        cls.def("getMapping", &PhotometryModel::getMapping, py::return_value_policy::reference_internal);
        cls.def("getTotalParameters", &PhotometryModel::getTotalParameters);

        cpputils::python::addOutputOp(cls, "__repr__");
        cls.def("__str__", [](PhotometryModel const &self) { return "PhotometryModel"; });
    });
}

void declareSimplePhotometryModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySimplePhotometryModel =
            py::class_<SimplePhotometryModel, std::shared_ptr<SimplePhotometryModel>, PhotometryModel>;

    wrappers.wrapType(PySimplePhotometryModel(wrappers.module, "SimplePhotometryModel"), [](auto &mod, auto &cls) {
        cls.def("__str__", [](SimplePhotometryModel const &self) { return "SimplePhotometryModel"; });
    });
}

void declareSimpleFluxModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySimpleFluxModel =
            py::class_<SimpleFluxModel, std::shared_ptr<SimpleFluxModel>, SimplePhotometryModel, PhotometryModel>;

    wrappers.wrapType(PySimpleFluxModel(wrappers.module, "SimpleFluxModel"), [](auto &mod, auto &cls) {
        cls.def(py::init<CcdImageList const &, double>(), "ccdImageList"_a, "errorPedestal"_a = 0);
        cls.def("__str__", [](SimpleFluxModel const &self) { return "SimpleFluxModel"; });
    });
}

void declareSimpleMagnitudeModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySimpleMagnitudeModel =  py::class_<SimpleMagnitudeModel, std::shared_ptr<SimpleMagnitudeModel>, SimplePhotometryModel,
               PhotometryModel>;

    wrappers.wrapType(PySimpleMagnitudeModel(wrappers.module, "SimpleMagnitudeModel"), [](auto &mod, auto &cls) {
        cls.def(py::init<CcdImageList const &, double>(), "ccdImageList"_a, "errorPedestal"_a = 0);
        cls.def("__str__", [](SimpleMagnitudeModel const &self) { return "SimpleMagnitudeModel"; });
    });
}

void declareConstrainedPhotometryModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyConstrainedPhotometryModel =
            py::class_<ConstrainedPhotometryModel, std::shared_ptr<ConstrainedPhotometryModel>, PhotometryModel>;

    wrappers.wrapType(
            PyConstrainedPhotometryModel(wrappers.module, "ConstrainedPhotometryModel"), [](auto &mod, auto &cls) {
                cls.def("__str__",
                        [](ConstrainedPhotometryModel const &self) { return "ConstrainedPhotometryModel"; });
            });
}

void declareConstrainedFluxModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyConstrainedFluxModel =  py::class_<ConstrainedFluxModel, std::shared_ptr<ConstrainedFluxModel>, PhotometryModel>;

    wrappers.wrapType(PyConstrainedFluxModel(wrappers.module, "ConstrainedFluxModel"), [](auto &mod, auto &cls) {
        cls.def(py::init<CcdImageList const &, lsst::geom::Box2D const &, int, double>(), "CcdImageList"_a,
                "bbox"_a, "visitOrder"_a = 7, "errorPedestal"_a = 0);
    });
}

void declareConstrainedMagnitudeModel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyConstrainedMagnitudeModel =  py::class_<ConstrainedMagnitudeModel, std::shared_ptr<ConstrainedMagnitudeModel>, PhotometryModel>;

    wrappers.wrapType(PyConstrainedMagnitudeModel(wrappers.module, "ConstrainedMagnitudeModel"), [](auto &mod, auto &cls) {
        cls.def(py::init<CcdImageList const &, lsst::geom::Box2D const &, int, double>(), "CcdImageList"_a,
                "bbox"_a, "visitOrder"_a = 7, "errorPedestal"_a = 0);
        cls.def("__str__", [](ConstrainedMagnitudeModel const &self) { return "ConstrainedMagnitudeModel"; });
    });
}
}  // namespace

void wrapPhotometryModels(lsst::cpputils::python::WrapperCollection &wrappers) {
    declarePhotometryModel(wrappers);
    declareSimplePhotometryModel(wrappers);
    declareSimpleFluxModel(wrappers);
    declareSimpleMagnitudeModel(wrappers);
    declareConstrainedPhotometryModel(wrappers);
    declareConstrainedFluxModel(wrappers);
    declareConstrainedMagnitudeModel(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
