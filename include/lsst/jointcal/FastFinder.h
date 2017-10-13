// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_FAST_FINDER_H
#define LSST_JOINTCAL_FAST_FINDER_H

#include <vector>
#include "lsst/jointcal/BaseStar.h"

namespace lsst {
namespace jointcal {

/*! \file
    \brief Fast locator in starlists.
*/

/*!  This is an auxillary class for matching objects from
  starlists. It allows to locate rapidly the closest objects from a
  given position. The very simple strategy is to sort objects
  according to 1 coordinate x, and to build an index that allows to
  select the objects with the x coordinate inside an interval. Then
  every slice in x is sorted according to y, which enables a fast scan
  inside a x slice. listMatchCollect takes about 10ms (PC 450 MHz,
  optimized "-O4") for a match between lists of about 2000 objects
  each, which is fast enough for our needs. The same "locator" is used
  in listMatchupShift, to avoid scanning the whole input lists. Timing
  on listMatchCollect and listMatchupShift indicates a gain in speed
  by more than one order of magnitude after implementation of this
  FastFinder.
*/

//! Fast locator in starlists.
class FastFinder {
public:
    const BaseStarList baselist;  // shallow copy of the initial list of stars (not used, acts as a
                                  // conservatory). The need is arguable.
    unsigned count;               // total number of objects (size of input list stars).
                                  /* the sorted pointer array: It does not seem very wise to use smart
                                     pointers here because reference counts will uselessly jump around
                                     during sorting */
    // Unfortunately we *must* use shared_ptr here because the return value of findClosest is used to pass
    // ownership
    std::vector<std::shared_ptr<const BaseStar>> stars;
    unsigned nslice;              // number of (X) slices
    std::vector<unsigned> index;  // index in "stars" of first object of each slice.
    double xmin, xmax, xstep;     // x bounds, slice size

    typedef decltype(stars)::value_type stars_element;
    typedef decltype(stars)::const_iterator pstar;

    //! Constructor
    FastFinder(const BaseStarList &list, const unsigned nXSlice = 100);

    //! Find the closest with some rejection capability
    std::shared_ptr<const BaseStar> findClosest(const Point &where, const double maxDist,
                                                bool (*SkipIt)(const BaseStar &) = nullptr) const;

    //!
    std::shared_ptr<const BaseStar> secondClosest(const Point &where, const double maxDist,
                                                  std::shared_ptr<const BaseStar> &closest,
                                                  bool (*SkipIt)(const BaseStar &) = nullptr) const;

    //! mostly for debugging
    void dump() const;

    //! Iterator meant to traverse objects within some limiting distance. Initializer is beginScan and end
    //! condition is (*it == NULL). Used by findClosest & co.

    class Iterator {
    public:  // could be made private, but what for??
        const FastFinder &finder;
        int currentSlice, endSlice;
        double yStart, yEnd;  // Y limits ( for all stripes)
        /* pointers to the first and beyond last stars in the y range for
           the current stripe :  */
        pstar current, pend;
        pstar null_value;  // with pointers being iterators, the null value is not NULL

        void check() const;

    public:
        Iterator(const FastFinder &f, const Point &where, double maxDist);
        void operator++();
        stars_element operator*() const;
    };

    Iterator beginScan(const Point &where, double maxDist) const;

    void findRangeInSlice(const int iSlice, const double yStart, const double yEnd, pstar &start,
                          pstar &end) const;
    pstar locateYStart(pstar begin, pstar end, double yVal) const;
    pstar locateYEnd(pstar begin, pstar end, double yVal) const;
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_FAST_FINDER_H
