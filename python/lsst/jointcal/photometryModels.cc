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

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotomModel.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/SimplePhotomModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotomModel(py::module &mod) {
    py::class_<PhotomModel, std::shared_ptr<PhotomModel>> cls(mod, "PhotomModel");

    cls.def("photomFactor", &SimplePhotomModel::photomFactor, "ccdImage"_a, "where"_a=Point());
}

void declareSimplePhotomModel(py::module &mod) {
    py::class_<SimplePhotomModel, std::shared_ptr<SimplePhotomModel>, PhotomModel> cls(mod,
                                                                                       "SimplePhotomModel");
    cls.def(py::init<CcdImageList const&>(), "ccdImageList"_a);
}

PYBIND11_PLUGIN(photometryModels) {
    py::module::import("lsst.jointcal.point"); // needed for photomFactor's default "where"
    py::module mod("photometryModels");

    declarePhotomModel(mod);
    declareSimplePhotomModel(mod);

    return mod.ptr();
}

}}}  // lsst::jointcal::<anonymous>
