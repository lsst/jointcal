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

#ifndef LSST_JOINTCAL_JOINTCAL_CONTROL_H
#define LSST_JOINTCAL_JOINTCAL_CONTROL_H

#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Source.h"

namespace lsst {
namespace jointcal {

struct JointcalControl {
    LSST_CONTROL_FIELD(sourceFluxField, std::string, "name of flux field in source catalog");

    explicit JointcalControl(std::string  sourceFluxField = "slot_CalibFlux")
            :  // Set sourceFluxType to the value used in the source selector.
              sourceFluxField(std::move(sourceFluxField)) {
        validate();
    }

    ~JointcalControl()= default;;

    void validate() const {
        if (sourceFluxField.empty()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "sourceFluxField must be specified");
        }
    }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_JOINTCAL_CONTROL_H
