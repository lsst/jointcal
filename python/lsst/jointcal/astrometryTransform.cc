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
#include "lsst/cpputils/python.h"

#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Point.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransform = py::class_<AstrometryTransform>;

    wrappers.wrapType(PyAstrometryTransform(wrappers.module, "AstrometryTransform"), [](auto &mod, auto &cls) {
        cls.def("apply",
                (jointcal::Point(AstrometryTransform::*)(const jointcal::Point &) const) &
                        AstrometryTransform::apply,
                "inPos"_a);
        cls.def("apply",
                (jointcal::Frame(AstrometryTransform::*)(Frame const &, bool) const) &AstrometryTransform::apply,
                "inputframe"_a, "inscribed"_a);
        cls.def("getNpar", &AstrometryTransform::getNpar);
        cls.def("offsetParams", &AstrometryTransform::offsetParams);
        cls.def("linearApproximation", &AstrometryTransform::linearApproximation);
        cls.def("computeDerivative", [](AstrometryTransform const &self, Point const &where, double const step) {
            AstrometryTransformLinear derivative;
            self.computeDerivative(where, derivative, step);
            return derivative;
        });

        cpputils::python::addOutputOp(cls, "__str__");
    });
}

void declareAstrometryTransformIdentity(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformIdentity =
            py::class_<AstrometryTransformIdentity, AstrometryTransform>;

    wrappers.wrapType(PyAstrometryTransformIdentity(wrappers.module, "AstrometryTransformIdentity"),
                      [](auto &mod, auto &cls) {});
}

void declareAstrometryTransformPolynomial(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformPolynomial =  py::class_<AstrometryTransformPolynomial,
               AstrometryTransform>;

    wrappers.wrapType(
            PyAstrometryTransformPolynomial(wrappers.module, "AstrometryTransformPolynomial"),
            [](auto &mod, auto &cls) {
                cls.def(py::init<const unsigned>(), "order"_a);
                cls.def("getOrder", &AstrometryTransformPolynomial::getOrder);
                cls.def("getCoefficient", py::overload_cast<std::size_t, std::size_t, std::size_t>(
                        &AstrometryTransformPolynomial::getCoefficient, py::const_));
                cls.def("setCoefficient",
                        [](AstrometryTransformPolynomial &self, std::size_t powX, std::size_t powY,
                           std::size_t whichCoord,
                           double value) { self.getCoefficient(powX, powY, whichCoord) = value; });
                cls.def("determinant", &AstrometryTransformPolynomial::determinant);
                cls.def("getNpar", &AstrometryTransformPolynomial::getNpar);
                cls.def("toAstMap", &AstrometryTransformPolynomial::toAstMap);
                cls.def("write", [](AstrometryTransformPolynomial const &self) {
                    std::stringstream result;
                    self.write(result);
                    return result.str();
                });
                cls.def("read", [](AstrometryTransformPolynomial &self, std::string const &str) {
                    std::istringstream istr(str);
                    self.read(istr);
                });
            });
}

void declareAstrometryTransformLinear(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformLinear = py::class_<AstrometryTransformLinear,
               AstrometryTransformPolynomial>;

    wrappers.wrapType(
            PyAstrometryTransformLinear(wrappers.module, "AstrometryTransformLinear"), [](auto &mod, auto &cls) {
                cls.def(py::init<>());
            });
}

void declareAstrometryTransformLinearShift(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformLinearShift = py::class_<AstrometryTransformLinearShift,
            AstrometryTransformLinear>;

    wrappers.wrapType(PyAstrometryTransformLinearShift(wrappers.module, "AstrometryTransformLinearShift"),
                      [](auto &mod, auto &cls) {});
}

void declareAstrometryTransformLinearRot(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformLinearRot =
            py::class_<AstrometryTransformLinearRot,
                    AstrometryTransformLinear>;

    wrappers.wrapType(
            PyAstrometryTransformLinearRot(wrappers.module, "AstrometryTransformLinearRot"),
            [](auto &mod, auto &cls) {});
}

void declareAstrometryTransformLinearScale(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformLinearScale =
            py::class_<AstrometryTransformLinearScale,
                    AstrometryTransformLinear>;

    wrappers.wrapType(
            PyAstrometryTransformLinearScale(wrappers.module, "AstrometryTransformLinearScale"),
            [](auto &mod, auto &cls) {});
}

void declareAstrometryTransformSkyWcs(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryTransformSkyWcs =
            py::class_<AstrometryTransformSkyWcs, AstrometryTransform>;

    wrappers.wrapType(
            PyAstrometryTransformSkyWcs(wrappers.module, "AstrometryTransformSkyWcs"), [](auto &mod, auto &cls) {
                cls.def("getSkyWcs", &AstrometryTransformSkyWcs::getSkyWcs);
            });
}

void declareBaseTanWcs(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyBaseTanWcs = py::class_<BaseTanWcs, AstrometryTransform> ;

    wrappers.wrapType(PyBaseTanWcs(wrappers.module, "BaseTanWcs"), [](auto &mod, auto &cls) {});
}

void declareTanPixelToRaDec(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyTanPixelToRaDec = py::class_<TanPixelToRaDec, AstrometryTransform>;

    wrappers.wrapType(PyTanPixelToRaDec(wrappers.module, "TanPixelToRaDec"), [](auto &mod, auto &cls) {});
}

void declareTanRaDecToPixel(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyTanRaDecToPixel = py::class_<TanRaDecToPixel, AstrometryTransform>;

    wrappers.wrapType(PyTanRaDecToPixel(wrappers.module, "TanRaDecToPixel"), [](auto &mod, auto &cls) {});
}

void declareTanSipPixelToRaDec(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyTanSipPixelToRaDec = py::class_<TanSipPixelToRaDec, BaseTanWcs>;

    wrappers.wrapType(PyTanSipPixelToRaDec(wrappers.module, "TanSipPixelToRaDec"), [](auto &mod, auto &cls) {});
}
}  // namespace

void wrapAstrometryTransform(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareAstrometryTransform(wrappers);
    declareAstrometryTransformIdentity(wrappers);
    declareAstrometryTransformPolynomial(wrappers);
    declareAstrometryTransformLinear(wrappers);
    declareAstrometryTransformLinearShift(wrappers);
    declareAstrometryTransformLinearRot(wrappers);
    declareAstrometryTransformLinearScale(wrappers);
    declareAstrometryTransformSkyWcs(wrappers);
    declareBaseTanWcs(wrappers);
    declareTanPixelToRaDec(wrappers);
    declareTanRaDecToPixel(wrappers);
    declareTanSipPixelToRaDec(wrappers);

    // utility functions
    wrappers.module.def("inversePolyTransform", &inversePolyTransform, "forward"_a, "domain"_a, "precision"_a,
            "maxOrder"_a = 9, "nSteps"_a = 50);
}

}  // namespace jointcal
}  // namespace lsst
