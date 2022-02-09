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

#ifndef LSST_JOINTCAL_TRIPLETLIST_H
#define LSST_JOINTCAL_TRIPLETLIST_H

#include "Eigen/Sparse"

#include <vector>

namespace lsst {
namespace jointcal {

using Trip = Eigen::Triplet<double>;

// at the moment this class implements the eigen format.
// it would be wise to implement it differently if talking to cholmod
class TripletList : public std::vector<Trip> {
public:
    TripletList(std::size_t count) : _nextFreeIndex(0) { reserve(count); };

    void addTriplet(Eigen::Index i, Eigen::Index j, double val) {
        push_back(Trip(i, j, val));
    }

    Eigen::Index getNextFreeIndex() const { return _nextFreeIndex; }

    void setNextFreeIndex(Eigen::Index index) { _nextFreeIndex = index; }

private:
    Eigen::Index _nextFreeIndex;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_TRIPLETLIST_H
