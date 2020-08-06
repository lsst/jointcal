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

#include "lsst/jointcal/CcdImage.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareCcdImage(py::module &mod) {
    py::class_<CcdImage, std::shared_ptr<CcdImage>> cls(mod, "CcdImage");

    cls.def(py::init<afw::table::SourceCatalog &, std::shared_ptr<lsst::afw::geom::SkyWcs>,
                     std::shared_ptr<lsst::afw::image::VisitInfo>, lsst::geom::Box2I const &,
                     std::string const &, std::shared_ptr<afw::image::PhotoCalib>,
                     std::shared_ptr<afw::cameraGeom::Detector>, int, int, std::string const &>(),
            "record"_a, "wcs"_a, "visitInfo"_a, "bbox"_a, "filter"_a, "photoCalib"_a, "detector"_a, "visit"_a,
            "ccd"_a, "fluxField"_a);

    cls.def("getPhotoCalib", &CcdImage::getPhotoCalib);

    cls.def("countStars", &CcdImage::countStars);

    cls.def("resetCatalogForFit", &CcdImage::resetCatalogForFit);

    cls.def("getBoresightRaDec", &CcdImage::getBoresightRaDec);
    cls.def_property_readonly("boresightRaDec", &CcdImage::getBoresightRaDec);

    cls.def("getCcdId", &CcdImage::getCcdId);
    cls.def_property_readonly("ccdId", &CcdImage::getCcdId);

    cls.def("getEpoch", &CcdImage::getEpoch);
    cls.def_property_readonly("epoch", &CcdImage::getEpoch);

    cls.def("getImageFrame", &CcdImage::getImageFrame, py::return_value_policy::reference_internal);
    cls.def_property_readonly("imageFrame", &CcdImage::getImageFrame,
                              py::return_value_policy::reference_internal);

    cls.def("getName", &CcdImage::getName);
    cls.def_property_readonly("name", &CcdImage::getName);

    cls.def("getVisit", &CcdImage::getVisit);
    cls.def_property_readonly("visit", &CcdImage::getVisit);

    cls.def("getDetector", &CcdImage::getDetector, py::return_value_policy::reference_internal);

    cls.def("getCommonTangentPoint", &CcdImage::getCommonTangentPoint,
            py::return_value_policy::reference_internal);
    cls.def("setCommonTangentPoint", &CcdImage::setCommonTangentPoint);
    cls.def_property("commonTangentPoint", &CcdImage::getCommonTangentPoint, &CcdImage::setCommonTangentPoint,
                     py::return_value_policy::reference_internal);

    cls.def("getSkyToTangentPlane", &CcdImage::getSkyToTangentPlane,
            py::return_value_policy::reference_internal);
    cls.def("getReadWcs", &CcdImage::getReadWcs, py::return_value_policy::reference_internal);
}

PYBIND11_MODULE(ccdImage, mod) { declareCcdImage(mod); }
}  // namespace
}  // namespace jointcal
}  // namespace lsst
