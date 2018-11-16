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

#include "lsst/jointcal/Gtransfo.h"
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

double StarMatch::computeChi2(const Gtransfo &gtransfo) const {
    FatPoint tr;
    gtransfo.transformPosAndErrors(point1, tr);
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

static std::unique_ptr<double[]> chi2_array(const StarMatchList &starMatchList, const Gtransfo &gtransfo) {
    unsigned s = starMatchList.size();
    auto res = std::unique_ptr<double[]>(new double[s]);
    unsigned count = 0;
    for (auto const &it : starMatchList) res[count++] = it.computeChi2(gtransfo);
    return res;
}

static unsigned chi2_cleanup(StarMatchList &starMatchList, const double chi2Cut, const Gtransfo &gtransfo) {
    unsigned erased = starMatchList.removeAmbiguities(gtransfo);
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
   If the transfo is not assigned, it will be set to a GtransfoLinear. User
   can set an other type/order using setTransfo() before call. */
void StarMatchList::refineTransfo(double nSigmas) {
    double cut;
    unsigned nremoved;
    if (!_transfo) _transfo.reset(new GtransfoLin);
    do {
        int nused = size();
        if (nused <= 2) {
            _chi2 = -1;
            break;
        }
        _chi2 = _transfo->fit(*this);
        /* convention of the fitted routines :
           -  chi2 = 0 means zero degrees of freedom
                (this was not enforced in Gtransfo{Lin,Quad,Cub} ...)
           -  chi2 = -1 means ndof <0 (and hence no possible fit)
           --> in either case, refinement is over
           The fact that chi2 = 0 was not enforced when necessary means
           that in this (rare) case, we were discarding matches at random....
           With GtransfoPoly::fit, this is no longer the case.
        */
        if (_chi2 <= 0) return;
        unsigned npair = int(size());
        if (npair == 0) break;  // should never happen

        // compute some chi2 statistics
        std::unique_ptr<double[]> chi2_array(new double[npair]);
        unsigned count = 0;
        for (auto &starMatch : *this) chi2_array[count++] = starMatch.chi2 = starMatch.computeChi2(*_transfo);

        std::sort(chi2_array.get(), chi2_array.get() + npair);
        double median = (npair & 1) ? chi2_array[npair / 2]
                                    : (chi2_array[npair / 2 - 1] + chi2_array[npair / 2]) * 0.5;

        // discard outliers : the cut is understood as a "distance" cut
        cut = sq(nSigmas) * median;
        nremoved = chi2_cleanup(*this, cut, *_transfo);
    } while (nremoved);
    _dist2 = computeDist2(*this, *_transfo);
}

/* not very robust : assumes that we went through Refine just before... */
double StarMatchList::computeResidual() const {
    int deno = (2. * size() - _transfo->getNpar());
    return (deno > 0) ? sqrt(_dist2 / deno) : -1;  // is -1 a good idea?
}

void StarMatchList::setDistance(const Gtransfo &gtransfo) {
    for (auto &smi : *this) smi.setDistance(gtransfo);  // c'est compact
}

unsigned StarMatchList::removeAmbiguities(const Gtransfo &gtransfo, int which) {
    if (!which) return 0;
    setDistance(gtransfo);
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

void StarMatchList::setTransfoOrder(int order) {
    if (order == 0)
        setTransfo(std::make_shared<GtransfoLinShift>());
    else if (order == 1)
        setTransfo(std::make_shared<GtransfoLin>());
    else
        setTransfo(GtransfoPoly(order));
    // might consider throwing if order does not make sense (e.g. >10)
    _order = order;
}

/* This routine should operate on a copy : refineTransfo
   might shorten the list */
/* it is not const although it tries not to change anything  */
std::unique_ptr<Gtransfo> StarMatchList::inverseTransfo() {
    if (!_transfo) return nullptr;

    auto old_transfo = _transfo->clone();
    double old_chi2 = _chi2;

    swap();
    setTransfoOrder(_order);
    refineTransfo(3.);  // keep same order
    auto inverted_transfo = _transfo->clone();
    setTransfo(old_transfo.get());
    swap();
    _chi2 = old_chi2;

    return inverted_transfo;
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
    GtransfoIdentity identity;
    for (auto const &starMatch : *this) {
        if (starMatch.computeDistance(identity) < mindist) n++;
    }
    return (n);
}

void StarMatchList::applyTransfo(StarMatchList &transformed, const Gtransfo *priorTransfo,
                                 const Gtransfo *posteriorTransfo) const {
    transformed.clear();
    GtransfoIdentity id;
    const Gtransfo &T1 = (priorTransfo) ? *priorTransfo : id;
    const Gtransfo &T2 = (posteriorTransfo) ? *posteriorTransfo : id;

    for (auto const &starMatch : *this) {
        FatPoint p1;
        T1.transformPosAndErrors(starMatch.point1, p1);
        FatPoint p2;
        T2.transformPosAndErrors(starMatch.point2, p2);
        transformed.push_back(StarMatch(p1, p2, starMatch.s1, starMatch.s2));
    }
}

void StarMatchList::dumpTransfo(std::ostream &stream) const {
    stream << " ================================================================" << std::endl
           << " Transformation between lists of order " << getTransfoOrder() << std::endl
           << *_transfo  //<< endl
           << " Chi2 = " << getChi2() << "  Residual = " << computeResidual() << std::endl
           << "  Number in the list = " << size() << std::endl
           << " ================================================================" << std::endl;
}

double computeDist2(const StarMatchList &starMatchList, const Gtransfo &gtransfo) {
    double dist2 = 0;
    for (auto const &starMatch : starMatchList)
        dist2 += gtransfo.apply(starMatch.point1).computeDist2(starMatch.point2);
    return dist2;
}

double computeChi2(const StarMatchList &starMatchList, const Gtransfo &gtransfo) {
    unsigned s = starMatchList.size();
    std::unique_ptr<double[]> chi2s(chi2_array(starMatchList, gtransfo));
    double chi2 = 0;
    for (unsigned k = 0; k < s; ++k) chi2 += chi2s[k];
    return chi2;
}
}  // namespace jointcal
}  // namespace lsst
