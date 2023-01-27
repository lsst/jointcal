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

namespace lsst {
namespace jointcal {

using lsst::cpputils::python::WrapperCollection;
void wrapAssociations(WrapperCollection &wrappers);
void wrapAstrometryMappings(WrapperCollection &wrappers);
void wrapAstrometryModels(WrapperCollection &wrappers);
void wrapAstrometryTransform(WrapperCollection &wrappers);
void wrapCcdImage(WrapperCollection &wrappers);
void wrapChi2(WrapperCollection &wrappers);
void wrapFitter(WrapperCollection &wrappers);
void wrapFrame(WrapperCollection &wrappers);
void wrapJointcalControl(WrapperCollection &wrappers);
void wrapPhotometryMappings(WrapperCollection &wrappers);
void wrapPhotometryModels(WrapperCollection &wrappers);
void wrapPhotometryTransform(WrapperCollection &wrappers);
void wrapProjectionHandler(WrapperCollection &wrappers);
void wrapStar(WrapperCollection &wrappers);

PYBIND11_MODULE(_jointcalLib, mod) {
    lsst::cpputils::python::WrapperCollection wrappers(mod, "lsst.jointcal");
    wrappers.addInheritanceDependency("lsst.sphgeom");
    wrapAssociations(wrappers);
    wrapAstrometryMappings(wrappers);
    wrapAstrometryModels(wrappers);
    wrapAstrometryTransform(wrappers);
    wrapCcdImage(wrappers);
    wrapChi2(wrappers);
    wrapFitter(wrappers);
    wrapFrame(wrappers);
    wrapJointcalControl(wrappers);
    wrapPhotometryMappings(wrappers);
    wrapPhotometryModels(wrappers);
    wrapPhotometryTransform(wrappers);
    wrapProjectionHandler(wrappers);
    wrapStar(wrappers);
    wrappers.finish();
}

}
}