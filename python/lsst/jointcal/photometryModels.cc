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
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/ConstrainedPhotometryModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryModel(py::module &mod) {
    py::class_<PhotometryModel, std::shared_ptr<PhotometryModel>> cls(mod, "PhotometryModel");

    cls.def("assignIndices", &PhotometryModel::assignIndices);
    cls.def("freezeErrorScales", &PhotometryModel::freezeErrorScales);
    cls.def("offsetParams", &PhotometryModel::offsetParams);
    cls.def("transform", &PhotometryModel::transform);
    cls.def("getMappingIndices", &PhotometryModel::getMappingIndices);
    cls.def("computeParameterDerivatives",
            [](PhotometryModel const &self, MeasuredStar const &star, CcdImage const &ccdImage) {
                Eigen::VectorXd derivatives(self.getNpar(ccdImage));
                self.computeParameterDerivatives(star, ccdImage, derivatives);
                return derivatives;
            });

    cls.def("getNpar", &PhotometryModel::getNpar);
    cls.def("toPhotoCalib", &PhotometryModel::toPhotoCalib);
    cls.def("__str__", [](PhotometryModel const &self) {
        std::ostringstream os;
        os << self;
        return os.str();
    });
}

void declareSimplePhotometryModel(py::module &mod) {
    py::class_<SimplePhotometryModel, std::shared_ptr<SimplePhotometryModel>, PhotometryModel> cls(
            mod, "SimplePhotometryModel");
    cls.def(py::init<CcdImageList const &>(), "ccdImageList"_a);
}

void declareConstrainedPhotometryModel(py::module &mod) {
    py::class_<ConstrainedPhotometryModel, std::shared_ptr<ConstrainedPhotometryModel>, PhotometryModel> cls(
            mod, "ConstrainedPhotometryModel");
    cls.def(py::init<CcdImageList const &, afw::geom::Box2D const &, int>(), "CcdImageList"_a, "bbox"_a,
            "visitDegree"_a = 7);
}

PYBIND11_PLUGIN(photometryModels) {
    py::module::import("lsst.jointcal.ccdImage");
    py::module::import("lsst.jointcal.photometryTransfo");
    py::module::import("lsst.jointcal.star");
    py::module mod("photometryModels");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declarePhotometryModel(mod);
    declareSimplePhotometryModel(mod);
    declareConstrainedPhotometryModel(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
