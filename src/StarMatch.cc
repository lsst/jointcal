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

#include <iostream>
#include <fstream>
#include <iomanip>

#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/BaseStar.h"
#include "algorithm"  // for copy

/* TO DO:
   think about imposing a maximum number of matches that may
   be discarded in Cleanup.
*/

namespace lsst {
namespace jointcal {

static double sq(double x) { return x * x; }

double StarMatch::computeChi2(const AstrometryTransform &transform) const {
    FatPoint tr;
    transform.transformPosAndErrors(point1, tr);
    double vxx = tr.vx + point2.vx;
    double vyy = tr.vy + point2.vy;
    double vxy = tr.vxy + point2.vxy;
    double det = vxx * vyy - vxy * vxy;
    return (vyy * sq(tr.x - point2.x) + vxx * sq(tr.y - point2.y) -
            2 * vxy * (tr.x - point2.x) * (tr.y - point2.y)) /
           det;
}

std::ostream &operator<<(std::ostream &stream, const StarMatch &match) {
    stream << match.point1.x << ' ' << match.point1.y << ' ' << match.point2.x << ' ' << match.point2.y << ' '
           << match.distance << std::endl;
    return stream;
}

std::ostream &operator<<(std::ostream &stream, const StarMatchList &starMatchList) {
    stream << " number of elements " << starMatchList.size() << std::endl;
    copy(starMatchList.begin(), starMatchList.end(), std::ostream_iterator<StarMatch>(stream));
    return stream;
}

static std::unique_ptr<double[]> chi2_array(const StarMatchList &starMatchList,
                                            const AstrometryTransform &transform) {
    unsigned s = starMatchList.size();
    auto res = std::unique_ptr<double[]>(new double[s]);
    unsigned count = 0;
    for (auto const &it : starMatchList) res[count++] = it.computeChi2(transform);
    return res;
}

static unsigned chi2_cleanup(StarMatchList &starMatchList, const double chi2Cut,
                             const AstrometryTransform &transform) {
    unsigned erased = starMatchList.removeAmbiguities(transform);
    for (auto smi = starMatchList.begin(); smi != starMatchList.end();) {
        if (smi->chi2 > chi2Cut) {
            smi = starMatchList.erase(smi);
            erased++;
        } else
            ++smi;
    }
    return erased;
}

/*! removes pairs beyond nSigmas in distance (where the sigma scale is
   set by the fit) and iterates until stabilization of the number of pairs.
   If the transform is not assigned, it will be set to a AstrometryTransformLinearear. User
   can set an other type/order using setTransform() before call. */
void StarMatchList::refineTransform(double nSigmas) {
    double cut;
    unsigned nremoved;
    if (!_transform) _transform.reset(new AstrometryTransformLinear);
    do {
        int nused = size();
        if (nused <= 2) {
            _chi2 = -1;
            break;
        }
        _chi2 = _transform->fit(*this);
        /* convention of the fitted routines :
           -  chi2 = 0 means zero degrees of freedom
                (this was not enforced in AstrometryTransform{Lin,Quad,Cub} ...)
           -  chi2 = -1 means ndof <0 (and hence no possible fit)
           --> in either case, refinement is over
           The fact that chi2 = 0 was not enforced when necessary means
           that in this (rare) case, we were discarding matches at random....
           With AstrometryTransformPolynomial::fit, this is no longer the case.
        */
        if (_chi2 <= 0) return;
        unsigned npair = int(size());
        if (npair == 0) break;  // should never happen

        // compute some chi2 statistics
        std::unique_ptr<double[]> chi2_array(new double[npair]);
        unsigned count = 0;
        for (auto &starMatch : *this)
            chi2_array[count++] = starMatch.chi2 = starMatch.computeChi2(*_transform);

        std::sort(chi2_array.get(), chi2_array.get() + npair);
        double median = (npair & 1) ? chi2_array[npair / 2]
                                    : (chi2_array[npair / 2 - 1] + chi2_array[npair / 2]) * 0.5;

        // discard outliers : the cut is understood as a "distance" cut
        cut = sq(nSigmas) * median;
        nremoved = chi2_cleanup(*this, cut, *_transform);
    } while (nremoved);
    _dist2 = computeDist2(*this, *_transform);
}

/* not very robust : assumes that we went through Refine just before... */
double StarMatchList::computeResidual() const {
    int deno = (2. * size() - _transform->getNpar());
    return (deno > 0) ? sqrt(_dist2 / deno) : -1;  // is -1 a good idea?
}

void StarMatchList::setDistance(const AstrometryTransform &transform) {
    for (auto &smi : *this) smi.setDistance(transform);  // c'est compact
}

unsigned StarMatchList::removeAmbiguities(const AstrometryTransform &transform, int which) {
    if (!which) return 0;
    setDistance(transform);
    int initial_count = size();
    if (which & 1) {
        sort(compareStar1);
        unique(sameStar1);
    }
    if (which & 2) {
        sort(compareStar2);
        unique(sameStar2);
    }
    return (initial_count - size());
}

void StarMatchList::setTransformOrder(int order) {
    if (order == 0)
        setTransform(std::make_shared<AstrometryTransformLinearShift>());
    else if (order == 1)
        setTransform(std::make_shared<AstrometryTransformLinear>());
    else
        setTransform(AstrometryTransformPolynomial(order));
    // might consider throwing if order does not make sense (e.g. >10)
    _order = order;
}

/* This routine should operate on a copy : refineTransform
   might shorten the list */
/* it is not const although it tries not to change anything  */
std::unique_ptr<AstrometryTransform> StarMatchList::inverseTransform() {
    if (!_transform) return nullptr;

    auto old_transform = _transform->clone();
    double old_chi2 = _chi2;

    swap();
    setTransformOrder(_order);
    refineTransform(3.);  // keep same order
    auto inverted_transform = _transform->clone();
    setTransform(old_transform.get());
    swap();
    _chi2 = old_chi2;

    return inverted_transform;
}

void StarMatchList::cutTail(int nKeep) {
    iterator si;
    int count = 0;
    for (si = begin(); si != end() && count < nKeep; ++count, ++si)
        ;
    erase(si, end());
}

void StarMatchList::swap() {
    for (auto &starMatch : *this) {
        starMatch.swap();
    }
}

int StarMatchList::recoveredNumber(double mindist) const {
    int n = 0;
    AstrometryTransformIdentity identity;
    for (auto const &starMatch : *this) {
        if (starMatch.computeDistance(identity) < mindist) n++;
    }
    return (n);
}

void StarMatchList::applyTransform(StarMatchList &transformed, const AstrometryTransform *priorTransform,
                                   const AstrometryTransform *posteriorTransform) const {
    transformed.clear();
    AstrometryTransformIdentity id;
    const AstrometryTransform &T1 = (priorTransform) ? *priorTransform : id;
    const AstrometryTransform &T2 = (posteriorTransform) ? *posteriorTransform : id;

    for (auto const &starMatch : *this) {
        FatPoint p1;
        T1.transformPosAndErrors(starMatch.point1, p1);
        FatPoint p2;
        T2.transformPosAndErrors(starMatch.point2, p2);
        transformed.push_back(StarMatch(p1, p2, starMatch.s1, starMatch.s2));
    }
}

void StarMatchList::printTransform(std::ostream &stream) const {
    stream << " ================================================================" << std::endl
           << " Transformation between lists of order " << getTransformOrder() << std::endl
           << *_transform  //<< endl
           << " Chi2 = " << getChi2() << "  Residual = " << computeResidual() << std::endl
           << "  Number in the list = " << size() << std::endl
           << " ================================================================" << std::endl;
}

double computeDist2(const StarMatchList &starMatchList, const AstrometryTransform &transform) {
    double dist2 = 0;
    for (auto const &starMatch : starMatchList)
        dist2 += transform.apply(starMatch.point1).computeDist2(starMatch.point2);
    return dist2;
}

double computeChi2(const StarMatchList &starMatchList, const AstrometryTransform &transform) {
    unsigned s = starMatchList.size();
    std::unique_ptr<double[]> chi2s(chi2_array(starMatchList, transform));
    double chi2 = 0;
    for (unsigned k = 0; k < s; ++k) chi2 += chi2s[k];
    return chi2;
}
}  // namespace jointcal
}  // namespace lsst
