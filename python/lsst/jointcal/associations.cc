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
#include "pybind11/stl.h"
#include "lsst/cpputils/python.h"

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/sphgeom/Circle.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAssociations(lsst::cpputils::python::WrapperCollection &wrappers) {
    using PyAssociations = py::class_<Associations>;
    wrappers.wrapType(PyAssociations(wrappers.module, "Associations"), [](auto &mod, auto &cls) {
        cls.def(py::init<>());
        cls.def(py::init<CcdImageList const &, double>(), "imageList"_a, "epoch"_a = 0);

        // NOTE: these could go away if the lists they wrap can be accessed directly.
        cls.def("refStarListSize", &Associations::refStarListSize);
        cls.def("fittedStarListSize", &Associations::fittedStarListSize);
        cls.def("associateCatalogs", &Associations::associateCatalogs, "matchCutInArcsec"_a = 0,
                "useFittedList"_a = false, "enlargeFittedList"_a = true);
        cls.def("collectRefStars", &Associations::collectRefStars, "refCat"_a, "matchCut"_a, "fluxField"_a,
                "refCoordinateErr"_a, "rejectBadFluxes"_a = false);
        cls.def("deprojectFittedStars", &Associations::deprojectFittedStars);
        cls.def("nCcdImagesValidForFit", &Associations::nCcdImagesValidForFit);
        cls.def("nFittedStarsWithAssociatedRefStar", &Associations::nFittedStarsWithAssociatedRefStar);

        cls.def("createCcdImage", &Associations::createCcdImage);
        cls.def("addCcdImage", &Associations::addCcdImage);
        cls.def("prepareFittedStars", &Associations::prepareFittedStars);
        cls.def("cleanFittedStars", &Associations::cleanFittedStars);

        cls.def("getCcdImageList", &Associations::getCcdImageList, py::return_value_policy::reference_internal);
        cls.def_property_readonly("ccdImageList", &Associations::getCcdImageList,
                                  py::return_value_policy::reference_internal);

        cls.def("computeBoundingCircle", &Associations::computeBoundingCircle);

        cls.def("getCommonTangentPoint", &Associations::getCommonTangentPoint);
        cls.def("setCommonTangentPoint", &Associations::setCommonTangentPoint);
        cls.def("computeCommonTangentPoint", &Associations::computeCommonTangentPoint);

        cls.def("getEpoch", &Associations::getEpoch);
        cls.def("setEpoch", &Associations::setEpoch);
    });
}
}  // namespace

void wrapAssociations(lsst::cpputils::python::WrapperCollection &wrappers) {
    declareAssociations(wrappers);
}
}  // namespace jointcal
}  // namespace lsst
