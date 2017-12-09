#include <algorithm>

#include "lsst/log/Log.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FastFinder.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.FastFinder");
}

namespace lsst {
namespace jointcal {

FastFinder::FastFinder(const BaseStarList &list, const unsigned nXSlice)
        : baselist(list), count(list.size()), stars(count), nslice(nXSlice), index(nslice + 1) {
    if (count == 0) return;

    // fill "stars"
    unsigned j = 0;
    for (auto const &ci : list) {
        stars[j] = ci;
        ++j;
    }

    sort(stars.begin(), stars.end(),
         [](const stars_element &E1, const stars_element &E2) { return (E1->x < E2->x); });

    xmin = stars[0]->x;
    xmax = stars[count - 1]->x;
    nslice = std::min(nslice, count);
    if (xmin == xmax) nslice = 1;

    // the x size of each slice:
    xstep = (xmax - xmin) / nslice;

    // fill the index array with the first star beyond the slice limit.
    index[0] = 0;  // first
    unsigned istar = 0;
    for (unsigned islice = 1; islice < nslice; ++islice) {
        double xend = xmin + (islice)*xstep;
        while (istar < count && stars[istar]->x < xend) ++istar;
        index[islice] = istar;
    }
    index[nslice] = count;  // last
    for (unsigned islice = 0; islice < nslice; ++islice) {
        sort(stars.begin() + index[islice], stars.begin() + index[islice + 1],
             [](const stars_element &E1, const stars_element &E2) {
                 return (E1->y < E2->y);
             });  // sort each slice in y.
    }
    // dump();
}

void FastFinder::dump() const {
    for (unsigned i = 0; i < count; ++i) {
        stars[i]->dump();
    }
}

std::shared_ptr<const BaseStar> FastFinder::findClosest(const Point &where, const double maxDist,
                                                        bool (*SkipIt)(const BaseStar &)) const {
    if (count == 0) return nullptr;
    FastFinder::Iterator it = beginScan(where, maxDist);
    if (*it == nullptr) return nullptr;
    std::shared_ptr<const BaseStar> pbest;
    double minDist2 = maxDist * maxDist;
    for (; *it != nullptr; ++it) {
        if (SkipIt && SkipIt(**it)) continue;
        double dist2 = where.computeDist2(**it);
        if (dist2 < minDist2) {
            pbest = *it;
            minDist2 = dist2;
        }
    }
    return pbest;
}

std::shared_ptr<const BaseStar> FastFinder::secondClosest(const Point &where, const double maxDist,
                                                          std::shared_ptr<const BaseStar> &closest,
                                                          bool (*SkipIt)(const BaseStar &)) const {
    closest = nullptr;
    if (count == 0) return nullptr;
    FastFinder::Iterator it = beginScan(where, maxDist);
    if (*it == nullptr) return nullptr;
    std::shared_ptr<const BaseStar> pbest1;  // closest
    std::shared_ptr<const BaseStar> pbest2;  // second closest
    double minDist1_2 = maxDist * maxDist;
    double minDist2_2 = maxDist * maxDist;
    for (; *it != nullptr; ++it) {
        if (SkipIt && SkipIt(**it)) continue;
        double dist2 = where.computeDist2(**it);
        if (dist2 < minDist1_2) {
            pbest2 = pbest1;
            minDist2_2 = minDist1_2;
            pbest1 = *it;
            minDist1_2 = dist2;
        } else if (dist2 < minDist2_2) {
            pbest2 = *it;
            minDist2_2 = dist2;
        }
    }
    closest = pbest1;
    return pbest2;
}

/* It is by no means clear the the 2 following routines are actually needed.
   It is nor clear to me (P.A) why they are different... but they really are.
*/
/* Locate the last position (in the sorted array) between begin and
   end that lies before yVal.*/
FastFinder::pstar FastFinder::locateYStart(pstar begin, pstar end, double yVal) const {
    if (begin == stars.end() || begin == end) return stars.end();
    int span = end - begin - 1;
    while (span > 1) {
        int half_span = span / 2;
        auto middle = begin + half_span;
        if ((*middle)->y < yVal) {
            begin += half_span;
            span -= half_span;
        } else {
            span -= (span - half_span);
        }
    }
    return begin;
}

/* Locate the first position (in the sorted array) between begin and
   end that lies beyond yVal.*/
FastFinder::pstar FastFinder::locateYEnd(pstar begin, pstar end, double yVal) const {
    if (begin == stars.end()) return stars.end();
    int span = end - begin - 1;
    while (span > 1) {
        int half_span = span / 2;
        auto middle = end - half_span;
        if ((*middle)->y > yVal) {
            end -= half_span;
            span -= half_span;
        } else {
            span -= (span - half_span);
        }
    }
    return end - 1;
}

void FastFinder::findRangeInSlice(const int iSlice, const double yStart, const double yEnd, pstar &start,
                                  pstar &end) const {
    start = locateYStart(stars.begin() + index[iSlice], stars.begin() + index[iSlice + 1], yStart);
    end = locateYEnd(start, stars.begin() + index[iSlice + 1], yEnd);
}

FastFinder::Iterator FastFinder::beginScan(const Point &where, double maxDist) const {
    return FastFinder::Iterator(*this, where, maxDist);
}

using Iterator = FastFinder::Iterator;

Iterator::Iterator(const FastFinder &F, const Point &where, double maxDist)
        : finder(F), null_value(F.stars.end()) {
    current = pend = null_value;  // does not iterate
    int startSlice = 0;
    if (finder.xstep != 0)  // means we have several slices
    {
        startSlice = std::max(0, int((where.x - maxDist - finder.xmin) / finder.xstep));
        /* obviously, endSlice (and starSlice) can be negative.
           This is why slice indices are "int" rather than "unsigned". */
        endSlice = std::min(int(finder.nslice), int((where.x + maxDist - finder.xmin) / finder.xstep) + 1);
    } else {
        startSlice = 0;
        endSlice = 1;
    }
    // beyond limits:
    if (startSlice >= int(finder.nslice) || endSlice < 0) return;
    // we are inside in x, so, we setup the y range:
    yStart = where.y - maxDist;
    yEnd = where.y + maxDist;
    /* rather than initializing here, we step back one
       slice and let "++" do its job */
    currentSlice = startSlice - 1;  // again, this requires "int" slices
    ++(*this);
}

FastFinder::stars_element Iterator::operator*() const {
    if (current != null_value) return *current;
    return nullptr;
}

void Iterator::operator++() {
    if (current != pend) {
        current++;
    } else {
        do {
            currentSlice++;
            if (currentSlice >= endSlice) {
                current = null_value;
                return;
            }
            finder.findRangeInSlice(currentSlice, yStart, yEnd, current, pend);
        } while (current == null_value);
    }
    check();
}

void FastFinder::Iterator::check() const {
    if (current != null_value &&
        (current < finder.stars.begin() || current >= finder.stars.begin() + finder.count)) {
        LOGLS_ERROR(_log, "Error in FastFinder " << *current << " " << *(finder.stars.begin()) << ' '
                                                 << *(finder.stars.begin() + finder.count));
    }
}
}  // namespace jointcal
}  // namespace lsst
