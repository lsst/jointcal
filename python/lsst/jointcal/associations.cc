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

    // NOTE: these could go away if the lists they wrap can be accessed directly.
    cls.def("refStarListSize", &Associations::refStarListSize);
    cls.def("fittedStarListSize", &Associations::fittedStarListSize);
    cls.def("associateCatalogs", &Associations::associateCatalogs, "matchCutInArcsec"_a = 0,
            "useFittedList"_a = false, "enlargeFittedList"_a = true);
    cls.def("collectRefStars", &Associations::collectRefStars);
    cls.def("deprojectFittedStars", &Associations::deprojectFittedStars);
    cls.def("nCcdImagesValidForFit", &Associations::nCcdImagesValidForFit);
    cls.def("nFittedStarsWithAssociatedRefStar", &Associations::nFittedStarsWithAssociatedRefStar);

    cls.def("addImage", &Associations::addImage);
    cls.def("selectFittedStars", &Associations::selectFittedStars);

    cls.def("getCcdImageList", &Associations::getCcdImageList, py::return_value_policy::reference_internal);
    cls.def_property_readonly("ccdImageList", &Associations::getCcdImageList,
                              py::return_value_policy::reference_internal);

    cls.def("getRaDecBBox", &Associations::getRaDecBBox);
    cls.def_property_readonly("raDecBBox", &Associations::getRaDecBBox);

    cls.def("getCommonTangentPoint", &Associations::getCommonTangentPoint);
    cls.def("setCommonTangentPoint", &Associations::setCommonTangentPoint);
    cls.def_property("commonTangentPoint", &Associations::getCommonTangentPoint,
                     &Associations::setCommonTangentPoint);
}

PYBIND11_PLUGIN(associations) {
    py::module::import("lsst.jointcal.ccdImage");
    py::module mod("associations");

    declareAssociations(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
