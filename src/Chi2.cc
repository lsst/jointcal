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

#include <utility>
#include <iostream>

#include "lsst/jointcal/Chi2.h"

namespace lsst {
namespace jointcal {

std::pair<double, double> Chi2List::computeAverageAndSigma() {
    double sum = 0;
    double sum2 = 0;
    for (const auto& i : *this) {
        sum += i.chi2;
        sum2 += std::pow(i.chi2, 2);
    }
    double average = sum / size();
    double sigma = sqrt(sum2 / size() - std::pow(average, 2));
    return std::make_pair(average, sigma);
}

std::ostream& operator<<(std::ostream& s, Chi2List const& chi2List) {
    s << "chi2 per star : ";
    for (const auto& chi2 : chi2List) {
        s << *(chi2.star) << " chi2: " << chi2.chi2 << " ; ";
    }
    s << std::endl;
    return s;
}

}  // namespace jointcal
}  // namespace lsst
