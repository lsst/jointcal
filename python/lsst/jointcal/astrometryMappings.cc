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

#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/ChipVisitAstrometryMapping.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryMapping(py::module &mod) {
    py::class_<AstrometryMapping, std::shared_ptr<AstrometryMapping>> cls(mod, "AstrometryMapping");

    cls.def("getNpar", &AstrometryMapping::getNpar);
    cls.def("transformPosAndErrors", [](AstrometryMapping const &self, jointcal::FatPoint &inPos) {
        jointcal::FatPoint outPos;
        self.transformPosAndErrors(inPos, outPos);
        return outPos;
    });
}

void declareChipVisitAstrometryMapping(py::module &mod) {
    py::class_<ChipVisitAstrometryMapping, std::shared_ptr<ChipVisitAstrometryMapping>, AstrometryMapping>
            cls(mod, "ChipVisitAstrometryMapping");

    cls.def("getTransform1", &ChipVisitAstrometryMapping::getTransform1,
            py::return_value_policy::reference_internal);
    cls.def_property_readonly("transform1", &ChipVisitAstrometryMapping::getTransform1,
                              py::return_value_policy::reference_internal);

    cls.def("getTransform2", &ChipVisitAstrometryMapping::getTransform2,
            py::return_value_policy::reference_internal);
    cls.def_property_readonly("transform2", &ChipVisitAstrometryMapping::getTransform2,
                              py::return_value_policy::reference_internal);
}

void declareSimpleAstrometryMapping(py::module &mod) {
    py::class_<SimpleAstrometryMapping, std::shared_ptr<SimpleAstrometryMapping>, AstrometryMapping> cls(
            mod, "SimpleAstrometryMapping");
    cls.def("getToBeFit", &SimpleAstrometryMapping::getToBeFit);
    cls.def("setToBeFit", &SimpleAstrometryMapping::setToBeFit);
    cls.def("getTransform", &SimpleAstrometryMapping::getTransform,
            py::return_value_policy::reference_internal);
}

void declareSimplePolyMapping(py::module &mod) {
    py::class_<SimplePolyMapping, std::shared_ptr<SimplePolyMapping>, SimpleAstrometryMapping> cls(
            mod, "SimplePolyMapping");
}

PYBIND11_MODULE(astrometryMappings, mod) {
    py::module::import("lsst.jointcal.star");
    py::module::import("lsst.jointcal.astrometryTransform");
    declareAstrometryMapping(mod);
    declareChipVisitAstrometryMapping(mod);
    declareSimpleAstrometryMapping(mod);
    declareSimplePolyMapping(mod);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
