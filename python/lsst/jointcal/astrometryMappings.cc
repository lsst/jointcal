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

#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/ChipVisitAstrometryMapping.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAstrometryMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAstrometryMapping = py::classh<AstrometryMapping>;
    wrappers.wrapType(PyAstrometryMapping(wrappers.module, "AstrometryMapping"), [](auto &mod, auto &cls) {
        cls.def("getNpar", &AstrometryMapping::getNpar);
        cls.def("transformPosAndErrors", [](AstrometryMapping const &self, jointcal::FatPoint &inPos) {
            jointcal::FatPoint outPos;
            self.transformPosAndErrors(inPos, outPos);
            return outPos;
        });
    });
}

void declareChipVisitAstrometryMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyChipVisitAstrometryMapping =
            py::classh<ChipVisitAstrometryMapping, AstrometryMapping>;

    wrappers.wrapType(PyChipVisitAstrometryMapping(wrappers.module, "ChipVisitAstrometryMapping"),
                      [](auto &mod, auto &cls) {
                          cls.def("getTransform1", &ChipVisitAstrometryMapping::getTransform1,
                                  py::return_value_policy::reference_internal);
                          cls.def_property_readonly("transform1", &ChipVisitAstrometryMapping::getTransform1,
                                                    py::return_value_policy::reference_internal);

                          cls.def("getTransform2", &ChipVisitAstrometryMapping::getTransform2,
                                  py::return_value_policy::reference_internal);
                          cls.def_property_readonly("transform2", &ChipVisitAstrometryMapping::getTransform2,
                                                    py::return_value_policy::reference_internal);
                      });
}

void declareSimpleAstrometryMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySimpleAstrometryMapping =
            py::classh<SimpleAstrometryMapping, AstrometryMapping>;
    wrappers.wrapType(PySimpleAstrometryMapping(wrappers.module, "impleAstrometryMapping"), [](auto &mod, auto &cls) {
        cls.def("getToBeFit", &SimpleAstrometryMapping::getToBeFit);
        cls.def("setToBeFit", &SimpleAstrometryMapping::setToBeFit);
        cls.def("getTransform", &SimpleAstrometryMapping::getTransform,
                py::return_value_policy::reference_internal);
    });
}

void declareSimplePolyMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PySimplePolyMapping = py::classh<SimplePolyMapping, SimpleAstrometryMapping>;

    wrappers.wrapType(PySimplePolyMapping(wrappers.module, "SimplePolyMapping"), [](auto &mod, auto &cls) {
    });
}
}  // namespace

void wrapAstrometryMappings(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareAstrometryMapping(wrappers);
    declareChipVisitAstrometryMapping(wrappers);
    declareSimpleAstrometryMapping(wrappers);
    declareSimplePolyMapping(wrappers);
}
}  // namespace jointcal
}  // namespace lsst
