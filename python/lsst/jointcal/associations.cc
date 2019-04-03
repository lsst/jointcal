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

#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/CcdImage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareAssociations(py::module &mod) {
    py::class_<Associations, std::shared_ptr<Associations>> cls(mod, "Associations");
    cls.def(py::init<>());
    cls.def(py::init<CcdImageList const &>(), "imageList"_a);

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
    cls.def("prepareFittedStars", &Associations::prepareFittedStars);

    cls.def("getCcdImageList", &Associations::getCcdImageList, py::return_value_policy::reference_internal);
    cls.def_property_readonly("ccdImageList", &Associations::getCcdImageList,
                              py::return_value_policy::reference_internal);

    cls.def("getRaDecBBox", &Associations::getRaDecBBox);
    cls.def_property_readonly("raDecBBox", &Associations::getRaDecBBox);

    cls.def("getCommonTangentPoint", &Associations::getCommonTangentPoint);
    cls.def("setCommonTangentPoint", &Associations::setCommonTangentPoint);
    cls.def("computeCommonTangentPoint", &Associations::computeCommonTangentPoint);
}

PYBIND11_MODULE(associations, mod) {
    py::module::import("lsst.jointcal.ccdImage");
    declareAssociations(mod);
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
