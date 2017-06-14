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

#include "lsst/jointcal/CcdImage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareCcdImage(py::module &mod) {
    py::class_<CcdImage, std::shared_ptr<CcdImage>> cls(mod, "CcdImage");

    // Not wrapping constructor: CcdImages are always built inside Associations.
    // cls.def(py::init<>());

    cls.def("getPhotoCalib", &CcdImage::getPhotoCalib);

    cls.def("getBoresightRaDec", &CcdImage::getBoresightRaDec);
    cls.def_property_readonly("boresightRaDec", &CcdImage::getBoresightRaDec);

    cls.def("getCcdId", &CcdImage::getCcdId);
    cls.def_property_readonly("ccdId", &CcdImage::getCcdId);

    cls.def("getImageFrame", &CcdImage::getImageFrame, py::return_value_policy::reference_internal);
    cls.def_property_readonly("imageFrame", &CcdImage::getImageFrame,
                              py::return_value_policy::reference_internal);

    cls.def("getVisit", &CcdImage::getVisit);
    cls.def_property_readonly("visit", &CcdImage::getVisit);

    cls.def("getCommonTangentPoint", &CcdImage::getCommonTangentPoint,
            py::return_value_policy::reference_internal);
    cls.def("setCommonTangentPoint", &CcdImage::setCommonTangentPoint);
    cls.def_property("commonTangentPoint", &CcdImage::getCommonTangentPoint, &CcdImage::setCommonTangentPoint,
                     py::return_value_policy::reference_internal);
}

PYBIND11_PLUGIN(ccdImage) {
    py::module mod("ccdImage");

    declareCcdImage(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
