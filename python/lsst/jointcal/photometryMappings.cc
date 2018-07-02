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
#include "pybind11/stl.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryMappingBase(py::module &mod) {
    py::class_<PhotometryMappingBase, std::shared_ptr<PhotometryMappingBase>> cls(mod,
                                                                                  "PhotometryMappingBase");

    cls.def("getNpar", &PhotometryMappingBase::getNpar);
    cls.def("transform", &PhotometryMappingBase::transform);
    cls.def("transformError", &PhotometryMappingBase::transformError);
    cls.def("computeParameterDerivatives",
            [](PhotometryMappingBase const &self, MeasuredStar const &star, double instFlux) {
                Eigen::VectorXd derivatives(self.getNpar());
                self.computeParameterDerivatives(star, instFlux, derivatives);
                return derivatives;
            });

    cls.def("offsetParams", &PhotometryMappingBase::offsetParams);
    cls.def("getParameters", &PhotometryMappingBase::getParameters);

    cls.def("getMappingIndices", [](PhotometryMappingBase const &self) {
        std::vector<unsigned> indices(0);
        self.getMappingIndices(indices);
        return indices;
    });

    cls.def("getIndex", &PhotometryMappingBase::getIndex);
    cls.def("setIndex", &PhotometryMappingBase::setIndex);
}

void declarePhotometryMapping(py::module &mod) {
    py::class_<PhotometryMapping, std::shared_ptr<PhotometryMapping>, PhotometryMappingBase> cls(
            mod, "PhotometryMapping");
    cls.def(py::init<std::shared_ptr<PhotometryTransfo>>(), "transfo"_a);

    cls.def("getTransfo", &PhotometryMapping::getTransfo);
    cls.def("getTransfoErrors", &PhotometryMapping::getTransfoErrors);
}

void declareChipVisitPhotometryMapping(py::module &mod) {
    py::class_<ChipVisitPhotometryMapping, std::shared_ptr<ChipVisitPhotometryMapping>, PhotometryMappingBase>
            cls(mod, "ChipVisitPhotometryMapping");

    cls.def(py::init<std::shared_ptr<PhotometryMapping>, std::shared_ptr<PhotometryMapping>, double>(),
            "chipMapping"_a, "visitMapping"_a, "err"_a);

    cls.def("getChipMapping", &ChipVisitPhotometryMapping::getChipMapping);
    cls.def("getVisitMapping", &ChipVisitPhotometryMapping::getVisitMapping);
}

PYBIND11_PLUGIN(photometryMappings) {
    py::module::import("lsst.jointcal.star");
    py::module::import("lsst.jointcal.photometryTransfo");
    py::module mod("photometryMappings");

    declarePhotometryMappingBase(mod);
    declarePhotometryMapping(mod);
    declareChipVisitPhotometryMapping(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
