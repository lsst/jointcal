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
#include "pybind11/eigen.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/utils/python.h"

#include "lsst/jointcal/PhotometryTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryTransform(py::module &mod) {
    py::class_<PhotometryTransform, std::shared_ptr<PhotometryTransform>> cls(mod, "PhotometryTransform");

    cls.def("transform",
            (double (PhotometryTransform::*)(double, double, double) const) & PhotometryTransform::transform,
            "x"_a, "y"_a, "value"_a);
    cls.def("transformError",
            (double (PhotometryTransform::*)(double, double, double, double) const) &
                    PhotometryTransform::transformError,
            "x"_a, "y"_a, "value"_a, "valueErr"_a);
    cls.def("offsetParams", &PhotometryTransform::offsetParams);
    cls.def("clone", &PhotometryTransform::clone);
    cls.def("getNpar", &PhotometryTransform::getNpar);
    cls.def("getParameters", &PhotometryTransform::getParameters);
    cls.def("computeParameterDerivatives",
            [](PhotometryTransform const &self, double x, double y, double instFlux) {
                Eigen::VectorXd derivatives(self.getNpar());
                self.computeParameterDerivatives(x, y, instFlux, derivatives);
                return derivatives;
            });

    utils::python::addOutputOp(cls, "__str__");
}

void declarePhotometryTransformSpatiallyInvariant(py::module &mod) {
    py::class_<PhotometryTransformSpatiallyInvariant, std::shared_ptr<PhotometryTransformSpatiallyInvariant>,
               PhotometryTransform>
            cls(mod, "PhotometryTransformSpatiallyInvariant");
}

void declareFluxTransfoSpatiallyInvariant(py::module &mod) {
    py::class_<FluxTransfoSpatiallyInvariant, std::shared_ptr<FluxTransfoSpatiallyInvariant>,
               PhotometryTransformSpatiallyInvariant, PhotometryTransform>
            cls(mod, "FluxTransfoSpatiallyInvariant");

    cls.def(py::init<double>(), "value"_a = 1);
}

void declareMagnitudeTransfoSpatiallyInvariant(py::module &mod) {
    py::class_<MagnitudeTransfoSpatiallyInvariant, std::shared_ptr<MagnitudeTransfoSpatiallyInvariant>,
               PhotometryTransformSpatiallyInvariant, PhotometryTransform>
            cls(mod, "MagnitudeTransfoSpatiallyInvariant");

    cls.def(py::init<double>(), "value"_a = 0);
}

void declarePhotometryTransformChebyshev(py::module &mod) {
    py::class_<PhotometryTransformChebyshev, std::shared_ptr<PhotometryTransformChebyshev>,
               PhotometryTransform>
            cls(mod, "PhotometryTransformChebyshev");

    cls.def("getCoefficients", &PhotometryTransformChebyshev::getCoefficients);
    cls.def("getOrder", &PhotometryTransformChebyshev::getOrder);
    cls.def("getBBox", &PhotometryTransformChebyshev::getBBox);
}

void declareFluxTransfoChebyshev(py::module &mod) {
    py::class_<FluxTransfoChebyshev, std::shared_ptr<FluxTransfoChebyshev>, PhotometryTransformChebyshev> cls(
            mod, "FluxTransfoChebyshev");

    cls.def(py::init<size_t, afw::geom::Box2D const &>(), "order"_a, "bbox"_a);
    cls.def(py::init<ndarray::Array<double, 2, 2> const &, afw::geom::Box2D const &>(), "coefficients"_a,
            "bbox"_a);
}

void declareMagnitudeTransfoChebyshev(py::module &mod) {
    py::class_<MagnitudeTransfoChebyshev, std::shared_ptr<MagnitudeTransfoChebyshev>,
               PhotometryTransformChebyshev>
            cls(mod, "MagnitudeTransfoChebyshev");

    cls.def(py::init<size_t, afw::geom::Box2D const &>(), "order"_a, "bbox"_a);
    cls.def(py::init<ndarray::Array<double, 2, 2> const &, afw::geom::Box2D const &>(), "coefficients"_a,
            "bbox"_a);
}

PYBIND11_MODULE(photometryTransform, mod) {
    declarePhotometryTransform(mod);

    declarePhotometryTransformSpatiallyInvariant(mod);
    declareFluxTransfoSpatiallyInvariant(mod);
    declareMagnitudeTransfoSpatiallyInvariant(mod);
    declarePhotometryTransformChebyshev(mod);
    declareFluxTransfoChebyshev(mod);
    declareMagnitudeTransfoChebyshev(mod);
}

}  // namespace
}  // namespace jointcal
}  // namespace lsst
