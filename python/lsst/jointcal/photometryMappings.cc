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

#include "lsst/jointcal/PhotometryTransform.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryMappingBase(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryMappingBase = py::classh<PhotometryMappingBase>;

    wrappers.wrapType(PyPhotometryMappingBase(wrappers.module, "PhotometryMappingBase"), [](auto &mod, auto &cls) {
        cls.def("getNpar", &PhotometryMappingBase::getNpar);
        cls.def("transform", &PhotometryMappingBase::transform);
        cls.def("transformError", &PhotometryMappingBase::transformError);
        cls.def("computeParameterDerivatives",
                [](PhotometryMappingBase const &self, MeasuredStar const &star, double instFlux) {
                    Eigen::VectorXd derivatives(self.getNpar());
                    self.computeParameterDerivatives(star, instFlux, derivatives);
                    return derivatives;
                });

        cls.def("getParameters", &PhotometryMappingBase::getParameters);

        cls.def("getMappingIndices", [](PhotometryMappingBase const &self) {
            IndexVector indices(0);
            self.getMappingIndices(indices);
            return indices;
        });

        cls.def("getIndex", &PhotometryMappingBase::getIndex);
        cls.def("setIndex", &PhotometryMappingBase::setIndex);
    });
}

void declarePhotometryMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyPhotometryMapping = py::classh<PhotometryMapping, PhotometryMappingBase>;

    wrappers.wrapType(PyPhotometryMapping(wrappers.module, "PhotometryMapping"), [](auto &mod, auto &cls) {
        cls.def(py::init<std::shared_ptr<PhotometryTransform>>(), "transform"_a);

        cls.def("offsetParams", &PhotometryMapping::offsetParams);
        cls.def("getTransform", &PhotometryMapping::getTransform);
        cls.def("getTransformErrors", &PhotometryMapping::getTransformErrors);
    });
}

void declareChipVisitPhotometryMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyChipVisitPhotometryMapping =
            py::classh<ChipVisitPhotometryMapping, PhotometryMappingBase>;

    wrappers.wrapType(
            PyChipVisitPhotometryMapping(wrappers.module, "ChipVisitPhotometryMapping"), [](auto &mod, auto &cls) {
                cls.def("setWhatToFit", &ChipVisitPhotometryMapping::setWhatToFit);
                cls.def("getChipMapping", &ChipVisitPhotometryMapping::getChipMapping);
                cls.def("getVisitMapping", &ChipVisitPhotometryMapping::getVisitMapping);
                cls.def("getNParChip", &ChipVisitPhotometryMapping::getNParChip);
                cls.def("getNParVisit", &ChipVisitPhotometryMapping::getNParVisit);
            });
}

void declareChipVisitFluxMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyChipVisitFluxMapping =
            py::classh<ChipVisitFluxMapping, ChipVisitPhotometryMapping>;

    wrappers.wrapType(PyChipVisitFluxMapping(wrappers.module, "ChipVisitFluxMapping"), [](auto &mod, auto &cls) {
        cls.def(py::init<std::shared_ptr<PhotometryMapping>, std::shared_ptr<PhotometryMapping>>(),
                "chipMapping"_a, "visitMapping"_a);
    });
}

void declareChipVisitMagnitudeMapping(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyChipVisitMagnitudeMapping = py::classh<ChipVisitMagnitudeMapping,
               ChipVisitPhotometryMapping>;

    wrappers.wrapType(PyChipVisitMagnitudeMapping(wrappers.module, "ChipVisitMagnitudeMapping"), [](auto &mod, auto &cls) {
        cls.def(py::init<std::shared_ptr<PhotometryMapping>, std::shared_ptr<PhotometryMapping>>(),
                "chipMapping"_a, "visitMapping"_a);
    });
}
}  // namespace

void wrapPhotometryMappings(lsst::cpputils::python::WrapperCollection &wrappers) {
    declarePhotometryMappingBase(wrappers);
    declarePhotometryMapping(wrappers);
    declareChipVisitPhotometryMapping(wrappers);
    declareChipVisitFluxMapping(wrappers);
    declareChipVisitMagnitudeMapping(wrappers);
}

}  // namespace jointcal
}  // namespace lsst
