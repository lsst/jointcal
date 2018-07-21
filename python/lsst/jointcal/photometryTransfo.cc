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
#include "pybind11/eigen.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/utils/python.h"

#include "lsst/jointcal/PhotometryTransfo.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryTransfo(py::module &mod) {
    py::class_<PhotometryTransfo, std::shared_ptr<PhotometryTransfo>> cls(mod, "PhotometryTransfo");

    cls.def("transform",
            (double (PhotometryTransfo::*)(double, double, double) const) & PhotometryTransfo::transform,
            "x"_a, "y"_a, "instFlux"_a);
    cls.def("offsetParams", &PhotometryTransfo::offsetParams);
    cls.def("clone", &PhotometryTransfo::clone);
    cls.def("getNpar", &PhotometryTransfo::getNpar);
    cls.def("getParameters", &PhotometryTransfo::getParameters);
    cls.def("computeParameterDerivatives",
            [](PhotometryTransfo const &self, double x, double y, double instFlux) {
                Eigen::VectorXd derivatives(self.getNpar());
                self.computeParameterDerivatives(x, y, instFlux, derivatives);
                return derivatives;
            });

    utils::python::addOutputOp(cls, "__str__");
}

void declarePhotometryTransfoSpatiallyInvariant(py::module &mod) {
    py::class_<PhotometryTransfoSpatiallyInvariant, std::shared_ptr<PhotometryTransfoSpatiallyInvariant>,
               PhotometryTransfo>
            cls(mod, "PhotometryTransfoSpatiallyInvariant");

    cls.def(py::init<double>(), "value"_a = 1);
}

void declarePhotometryTransfoChebyshev(py::module &mod) {
    py::class_<PhotometryTransfoChebyshev, std::shared_ptr<PhotometryTransfoChebyshev>, PhotometryTransfo>
            cls(mod, "PhotometryTransfoChebyshev");

    cls.def(py::init<size_t, afw::geom::Box2D const &>(), "order"_a, "bbox"_a);
    cls.def(py::init<ndarray::Array<double, 2, 2> const &, afw::geom::Box2D const &>(), "coefficients"_a,
            "bbox"_a);

    cls.def("getCoefficients", &PhotometryTransfoChebyshev::getCoefficients);
    cls.def("getOrder", &PhotometryTransfoChebyshev::getOrder);
    cls.def("getBBox", &PhotometryTransfoChebyshev::getBBox);
}

PYBIND11_MODULE(photometryTransfo, mod) {
    declarePhotometryTransfo(mod);
    declarePhotometryTransfoSpatiallyInvariant(mod);
    declarePhotometryTransfoChebyshev(mod);
}

}  // namespace
}  // namespace jointcal
}  // namespace lsst
