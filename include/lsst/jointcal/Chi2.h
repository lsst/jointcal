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

#ifndef LSST_JOINTCAL_CHI2_H
#define LSST_JOINTCAL_CHI2_H

#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "lsst/jointcal/BaseStar.h"

namespace lsst {
namespace jointcal {

/**
 * Base class for Chi2Statistic and Chi2List, to allow addEntry inside Fitter for either class.
 *
 * Essentially a mixin.
 */
class Chi2Accumulator {
public:
    virtual void addEntry(double inc, std::ptrdiff_t dof, std::shared_ptr<BaseStar> star) = 0;

    virtual ~Chi2Accumulator(){};
};

/// Simple structure to accumulate chi2 and ndof.
class Chi2Statistic : public Chi2Accumulator {
public:
    double chi2;
    std::ptrdiff_t ndof;

    Chi2Statistic() : chi2(0), ndof(0){};

    friend std::ostream& operator<<(std::ostream& s, Chi2Statistic const& chi2) {
        s << "chi2/ndof : " << chi2.chi2 << '/' << chi2.ndof << '=' << chi2.chi2 / chi2.ndof;
        return s;
    }

    // Addentry has an ignored third argument in order to make it compatible with Chi2List.
    void addEntry(double inc, std::ptrdiff_t dof, std::shared_ptr<BaseStar>) override {
        chi2 += inc;
        ndof += dof;
    }

    Chi2Statistic& operator+=(Chi2Statistic const& rhs) {
        chi2 += rhs.chi2;
        ndof += rhs.ndof;
        return *this;
    }
};

/*
 * A class to accumulate chi2 contributions together with pointers to the contributors.
 *
 * This structure lets one compute the chi2 statistics (average and variance) and directly point back
 * to the bad guys without relooping.
 * The Chi2Star routine makes it compatible with AstrometryFit's
 * accumulateStatImage and accumulateStatImageList.
 */
struct Chi2Star {
    double chi2;
    std::shared_ptr<BaseStar> star;

    Chi2Star(double chi2, std::shared_ptr<BaseStar> star) : chi2(chi2), star(std::move(star)) {}
    // for sorting
    bool operator<(Chi2Star const& rhs) const { return (chi2 < rhs.chi2); }

    friend std::ostream& operator<<(std::ostream& s, Chi2Star const& chi2Star) {
        s << "chi2: " << chi2Star.chi2 << " star: " << *(chi2Star.star) << std::endl;
        return s;
    }
};

/// Structure to accumulate the chi2 contributions per each star (to help find outliers).
class Chi2List : public Chi2Accumulator, public std::vector<Chi2Star> {
public:
    void addEntry(double chi2, std::ptrdiff_t ndof, std::shared_ptr<BaseStar> star) override {
        push_back(Chi2Star(chi2, std::move(star)));
    }

    /// Compute the average and std-deviation of these chisq values.
    std::pair<double, double> computeAverageAndSigma();

    friend std::ostream& operator<<(std::ostream& s, Chi2List const& chi2List);
};

}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_CHI2_H
