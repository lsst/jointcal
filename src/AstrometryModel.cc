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

#include "lsst/jointcal/AstrometryModel.h"

namespace lsst {
namespace jointcal {

bool AstrometryModel::validate(CcdImageList const &ccdImageList, int ndof) const {
    bool check = true;
    if (ndof < 0) {
        check &= false;
        LOGLS_ERROR(_log, "This model only has "
                                  << ndof << " degrees of freedom, with " << getTotalParameters()
                                  << " total parameters. Reduce the model complexity (e.g. polynomial order)"
                                     " to better match the number of measured sources.");
    }
    return check;
}

std::ostream &operator<<(std::ostream &stream, AstrometryModel const &model) {
    stream << model.toString();
    return stream;
}

}  // namespace jointcal
}  // namespace lsst
