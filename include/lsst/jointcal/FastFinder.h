#ifndef FASTFINDER__H
#define FASTFINDER__H

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
  inside a x slice. ListMatchCollect takes about 10ms (PC 450 MHz,
  optimized "-O4") for a match between lists of about 2000 objects
  each, which is fast enough for our needs. The same "locator" is used
  in ListMatchupShift, to avoid scanning the whole input lists. Timing
  on ListMatchCollect and ListMatchupShift indicates a gain in speed
  by more than one order of magnitude after implementation of this
  FastFinder.
*/

//! Fast locator in starlists.
class FastFinder
{
  //  private :
  public :
  const BaseStarList baselist; // shallow copy of the initial list of stars (not used, acts as a conservatory). The need is arguable.
  unsigned count;  // total number of objects (size of input list stars).
  /* the sorted pointer array: It does not seem very wise to use smart
     pointers here because reference counts will uselessly jump around
     during sorting */
  std::vector<const BaseStar*>  stars;
  unsigned nslice; // number of (X) slices
  std::vector<unsigned> index;// index in "stars" of first object of each slice.
  double xmin,xmax, xstep; // x bounds, slice size

  typedef decltype(stars)::value_type stars_element;
  typedef decltype(stars)::const_iterator pstar;

public :
    //! Constructor
  FastFinder(const BaseStarList &List, const unsigned NXslice = 100);

  //! Find the closest with some rejection capability
  const BaseStar *FindClosest(const Point &Where, const double MaxDist,
			      bool (*SkipIt)(const BaseStar *) = NULL) const;

  //!
  const BaseStar *SecondClosest(const Point &Where,
				const double MaxDist,
				const BaseStar* &Closest,
  				bool (*SkipIt)(const BaseStar *)= NULL) const;


  //! mostly for debugging
  void dump() const;


  //! Iterator meant to traverse objects within some limiting distance. Initializer is begin_scan and end condition is (*it == NULL). Used by FindClosest & co.

  class Iterator {
    public : // could be made private, but what for??
    const FastFinder &finder;
    int currentSlice, endSlice;
    double yStart, yEnd; // Y limits ( for all stripes)
    /* pointers to the first and beyond last stars in the y range for
       the current stripe :  */
    pstar current, pend;
    pstar null_value; // with pointers being iterators, the null value is not NULL

    void check() const;

  public:
    Iterator(const FastFinder &f, const Point &Where, double MaxDist);
    void operator++() ;
    stars_element operator*() const;
  };

#ifndef SWIG
  Iterator begin_scan(const Point &Where, double MaxDist) const;
#endif

  void find_range_in_slice(const int iSlice, const double YStart, const double YEnd, pstar &Start, pstar &End) const;
  pstar locate_y_start(pstar Begin, pstar End, double YVal) const;
  pstar locate_y_end(pstar Begin, pstar End, double YVal) const;

};

}}
#endif /* FASTFINDER__H */
