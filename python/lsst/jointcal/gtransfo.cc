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

#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/SipToGtransfo.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declareGtransfo(py::module &mod) {
    py::class_<Gtransfo, std::shared_ptr<Gtransfo>> cls(mod, "Gtransfo");
}

void declareBaseTanWcs(py::module &mod) {
    py::class_<BaseTanWcs, std::shared_ptr<BaseTanWcs>, Gtransfo> cls(mod, "BaseTanWcs");
}

void declareTanRaDec2Pix(py::module &mod) {
    py::class_<TanRaDec2Pix, std::shared_ptr<TanRaDec2Pix>, Gtransfo> cls(mod, "TanRaDec2Pix");
}

void declareTanSipPix2RaDec(py::module &mod) {
    py::class_<TanSipPix2RaDec, std::shared_ptr<TanSipPix2RaDec>, BaseTanWcs> cls(mod, "TanSipPix2RaDec");
}

PYBIND11_PLUGIN(gtransfo) {
    py::module::import("lsst.jointcal.frame");
    py::module mod("gtransfo");

    declareGtransfo(mod);
    declareBaseTanWcs(mod);
    declareTanRaDec2Pix(mod);
    declareTanSipPix2RaDec(mod);

    // utility functions
    mod.def("gtransfoToTanWcs", &gtransfoToTanWcs); // from SipToGtransfo.h

    return mod.ptr();
}

}}}  // lsst::jointcal::<anonymous>
