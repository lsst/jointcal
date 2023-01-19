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

#include "astshim.h"
#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/utils/python.h"

#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Point.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryTransform(py::module &mod) {
    py::class_<AstrometryTransform, std::shared_ptr<AstrometryTransform>> cls(mod, "AstrometryTransform");

    cls.def("apply",
            (jointcal::Point(AstrometryTransform::*)(const jointcal::Point &) const) &
                    AstrometryTransform::apply,
            "inPos"_a);
    cls.def("apply",
            (jointcal::Frame(AstrometryTransform::*)(Frame const &, bool) const) & AstrometryTransform::apply,
            "inputframe"_a, "inscribed"_a);
    cls.def("getNpar", &AstrometryTransform::getNpar);
    cls.def("offsetParams", &AstrometryTransform::offsetParams);
    cls.def("linearApproximation", &AstrometryTransform::linearApproximation);
    cls.def("computeDerivative", [](AstrometryTransform const &self, Point const &where, double const step) {
        AstrometryTransformLinear derivative;
        self.computeDerivative(where, derivative, step);
        return derivative;
    });

    utils::python::addOutputOp(cls, "__str__");
}

void declareAstrometryTransformIdentity(py::module &mod) {
    py::class_<AstrometryTransformIdentity, std::shared_ptr<AstrometryTransformIdentity>, AstrometryTransform>
            cls(mod, "AstrometryTransformIdentity");
}

void declareAstrometryTransformPolynomial(py::module &mod) {
    py::class_<AstrometryTransformPolynomial, std::shared_ptr<AstrometryTransformPolynomial>,
               AstrometryTransform>
            cls(mod, "AstrometryTransformPolynomial");

    cls.def(py::init<const unsigned>(), "order"_a);
    cls.def("getOrder", &AstrometryTransformPolynomial::getOrder);
    cls.def("getCoefficient", py::overload_cast<std::size_t, std::size_t, std::size_t>(
                                      &AstrometryTransformPolynomial::getCoefficient, py::const_));
    cls.def("setCoefficient", [](AstrometryTransformPolynomial &self, std::size_t powX, std::size_t powY,
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
}

void declareAstrometryTransformLinear(py::module &mod) {
    py::class_<AstrometryTransformLinear, std::shared_ptr<AstrometryTransformLinear>,
               AstrometryTransformPolynomial>
            cls(mod, "AstrometryTransformLinear");

    cls.def(py::init<>());
}

void declareAstrometryTransformLinearShift(py::module &mod) {
    py::class_<AstrometryTransformLinearShift, std::shared_ptr<AstrometryTransformLinearShift>,
               AstrometryTransformLinear>
            cls(mod, "AstrometryTransformLinearShift");
}

void declareAstrometryTransformLinearRot(py::module &mod) {
    py::class_<AstrometryTransformLinearRot, std::shared_ptr<AstrometryTransformLinearRot>,
               AstrometryTransformLinear>
            cls(mod, "AstrometryTransformLinearRot");
}

void declareAstrometryTransformLinearScale(py::module &mod) {
    py::class_<AstrometryTransformLinearScale, std::shared_ptr<AstrometryTransformLinearScale>,
               AstrometryTransformLinear>
            cls(mod, "AstrometryTransformLinearScale");
}

void declareAstrometryTransformSkyWcs(py::module &mod) {
    py::class_<AstrometryTransformSkyWcs, std::shared_ptr<AstrometryTransformSkyWcs>, AstrometryTransform>
            cls(mod, "AstrometryTransformSkyWcs");
    cls.def("getSkyWcs", &AstrometryTransformSkyWcs::getSkyWcs);
}

void declareBaseTanWcs(py::module &mod) {
    py::class_<BaseTanWcs, std::shared_ptr<BaseTanWcs>, AstrometryTransform> cls(mod, "BaseTanWcs");
}

void declareTanPixelToRaDec(py::module &mod) {
    py::class_<TanPixelToRaDec, std::shared_ptr<TanPixelToRaDec>, AstrometryTransform> cls(mod,
                                                                                           "TanPixelToRaDec");
}

void declareTanRaDecToPixel(py::module &mod) {
    py::class_<TanRaDecToPixel, std::shared_ptr<TanRaDecToPixel>, AstrometryTransform> cls(mod,
                                                                                           "TanRaDecToPixel");
}

void declareTanSipPixelToRaDec(py::module &mod) {
    py::class_<TanSipPixelToRaDec, std::shared_ptr<TanSipPixelToRaDec>, BaseTanWcs> cls(mod,
                                                                                        "TanSipPixelToRaDec");
}

PYBIND11_MODULE(astrometryTransform, mod) {
    py::module::import("astshim");
    py::module::import("lsst.jointcal.frame");
    py::module::import("lsst.jointcal.star");
    declareAstrometryTransform(mod);
    declareAstrometryTransformIdentity(mod);
    declareAstrometryTransformPolynomial(mod);
    declareAstrometryTransformLinear(mod);
    declareAstrometryTransformLinearShift(mod);
    declareAstrometryTransformLinearRot(mod);
    declareAstrometryTransformLinearScale(mod);
    declareAstrometryTransformSkyWcs(mod);
    declareBaseTanWcs(mod);
    declareTanPixelToRaDec(mod);
    declareTanRaDecToPixel(mod);
    declareTanSipPixelToRaDec(mod);

    // utility functions
    mod.def("inversePolyTransform", &inversePolyTransform, "forward"_a, "domain"_a, "precision"_a,
            "maxOrder"_a = 9, "nSteps"_a = 50);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
