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
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/jointcal/PhotometryTransform.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryTransform = py::class_<PhotometryTransform, std::shared_ptr<PhotometryTransform>>;

    wrappers.wrapType(PyPhotometryTransform(wrappers.module, "PhotometryTransform"), [](auto &mod, auto &cls) {
        cls.def("transform",
                (double (PhotometryTransform::*)(double, double, double) const) &PhotometryTransform::transform,
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
        cpputils::python::addOutputOp(cls, "__str__");
        cpputils::python::addOutputOp(cls, "__repr__");
    });
}

void declarePhotometryTransformSpatiallyInvariant(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryTransformSpatiallyInvariant = py::class_<PhotometryTransformSpatiallyInvariant, std::shared_ptr<PhotometryTransformSpatiallyInvariant>,
               PhotometryTransform>;

    wrappers.wrapType(
            PyPhotometryTransformSpatiallyInvariant(
                    wrappers.module, "PhotometryTransformSpatiallyInvariant"), [](auto &mod, auto &cls) {});
}

void declareFluxTransformSpatiallyInvariant(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyFluxTransformSpatiallyInvariant = py::class_<FluxTransformSpatiallyInvariant, std::shared_ptr<FluxTransformSpatiallyInvariant>,
               PhotometryTransformSpatiallyInvariant, PhotometryTransform>;

    wrappers.wrapType(
            PyFluxTransformSpatiallyInvariant(
                    wrappers.module, "FluxTransformSpatiallyInvariant"), [](auto &mod, auto &cls) {
                cls.def(py::init<double>(), "value"_a = 1);
            });
}

void declareMagnitudeTransformSpatiallyInvariant(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyMagnitudeTransformSpatiallyInvariant =  py::class_<MagnitudeTransformSpatiallyInvariant, std::shared_ptr<MagnitudeTransformSpatiallyInvariant>,
               PhotometryTransformSpatiallyInvariant, PhotometryTransform>;

    wrappers.wrapType(
            PyMagnitudeTransformSpatiallyInvariant(
                    wrappers.module, "MagnitudeTransformSpatiallyInvariant"), [](auto &mod, auto &cls) {
                cls.def(py::init<double>(), "value"_a = 0);
            });
}

void declarePhotometryTransformChebyshev(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryTransformChebyshev =
            py::class_<PhotometryTransformChebyshev, std::shared_ptr<PhotometryTransformChebyshev>, PhotometryTransform>;

    wrappers.wrapType(PyPhotometryTransformChebyshev(
            wrappers.module, "PhotometryTransformChebyshev"), [](auto &mod, auto &cls) {
        cls.def("getCoefficients", &PhotometryTransformChebyshev::getCoefficients);
        cls.def("getOrder", &PhotometryTransformChebyshev::getOrder);
        cls.def("getBBox", &PhotometryTransformChebyshev::getBBox);
        cls.def("integrate", py::overload_cast<>(&PhotometryTransformChebyshev::integrate, py::const_));
        cls.def("integrate",
                py::overload_cast<geom::Box2D const &>(&PhotometryTransformChebyshev::integrate, py::const_),
                "box"_a);
    });
}

void declareFluxTransformChebyshev(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyFluxTransformChebyshev =
            py::class_<FluxTransformChebyshev, std::shared_ptr<FluxTransformChebyshev>, PhotometryTransformChebyshev>;

    wrappers.wrapType(PyFluxTransformChebyshev(wrappers.module, "FluxTransformChebyshev"), [](auto &mod, auto &cls) {
        cls.def(py::init<size_t, geom::Box2D const &>(), "order"_a, "bbox"_a);
        cls.def(py::init<ndarray::Array<double, 2, 2> const &, geom::Box2D const &>(), "coefficients"_a,
                "bbox"_a);
    });
}

void declareMagnitudeTransformChebyshev(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyMagnitudeTransformChebyshev = py::class_<MagnitudeTransformChebyshev, std::shared_ptr<MagnitudeTransformChebyshev>,
               PhotometryTransformChebyshev>;

    wrappers.wrapType(
            PyMagnitudeTransformChebyshev(wrappers.module, "MagnitudeTransformChebyshev"), [](auto &mod, auto &cls) {
                cls.def(py::init<size_t, geom::Box2D const &>(), "order"_a, "bbox"_a);
                cls.def(py::init<ndarray::Array<double, 2, 2> const &, geom::Box2D const &>(), "coefficients"_a,
                        "bbox"_a);
            });
}
}  // namespace

void wrapPhotometryTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    declarePhotometryTransform(wrappers);
    declarePhotometryTransformSpatiallyInvariant(wrappers);
    declareFluxTransformSpatiallyInvariant(wrappers);
    declareMagnitudeTransformSpatiallyInvariant(wrappers);
    declarePhotometryTransformChebyshev(wrappers);
    declareFluxTransformChebyshev(wrappers);
    declareMagnitudeTransformChebyshev(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
