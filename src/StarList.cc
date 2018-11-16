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

#ifndef STARLIST__CC
#define STARLIST__CC

#include "lsst/pex/exceptions.h"

#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

template <class Star>
void StarList<Star>::fluxSort() {
    typedef StarList<Star>::Element E;
    this->sort([](const E &e1, const E &e2) { return (e1->getFlux() > e2->getFlux()); });
}

template <class Star>
void StarList<Star>::cutTail(const int nKeep) {
    int count = 0;
    auto si = this->begin();
    for (; si != this->end() && count < nKeep; ++count, ++si)
        ;
    while (si != this->end()) {
        si = this->erase(si);
    }
}

template <class Star>
void StarList<Star>::extractInFrame(StarList<Star> &out, const Frame &frame) const {
    for (auto const &star : *this) {
        if (frame.inFrame(*star)) {
            out.push_back(std::make_shared<Star>(*star));
        }
    }
}

template <class Star>
void StarList<Star>::copyTo(StarList<Star> &copy) const {
    copy.clearList();
    for (auto const &si : *this) copy.push_back(std::make_shared<Star>(*si));
}

// Explicit instantiations
template class StarList<BaseStar>;
template class StarList<FittedStar>;
template class StarList<MeasuredStar>;
}  // namespace jointcal
}  // namespace lsst

#endif /* STARLIST__CC */
