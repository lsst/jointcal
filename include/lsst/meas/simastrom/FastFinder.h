#ifndef FASTFINDER__H
#define FASTFINDER__H

#include <vector>
#include "lsst/meas/simastrom/BaseStar.h"

namespace lsst {
namespace meas {
namespace simastrom {

/*! \file
    \brief Fast locator in starlists.
*/
 
    /*!
  This is an auxillary class for  matching objects from starlists. It allows to locate rapidly
  the closest objects from a given position. The very simple strategy is to sort objects according
  to 1 coordinate x, and to build an index that allows to select the objects with the
  x coordinate inside an interval. Then every slice in x is sorted according to y, which enables 
  a fast scan inside a x slice. ListMatchCollect takes about 10ms (PC 450 MHz, optimized "-O4") 
  for a  match between lists of about 2000 objects each, which is fast enough for our needs.
  The same "locator" is used in ListMatchupShift, to avoid scanning the whole input lists.
  Timing on ListMatchCollect and ListMatchupShift indicates a gain in speed by more than one order
  after implementation of this FastFinder. One should not delete objects in the list passed to the 
  FastFinder constructor
  during the whole life of the FastFinder object.
*/

//! Fast locator in starlists.
class FastFinder
{
  //  private :
  public :
  const BaseStarList baselist; // pointer to the initial list of stars
  unsigned count;               // total number of objects (size of stars).
  std::vector<const BaseStar*>  stars;
  typedef decltype(stars)::value_type stars_element;
  typedef decltype(stars)::const_iterator pstar;
  std::vector<unsigned> index;              // index in "stars" of first object of each slice.
  double xmin,xmax, xstep; // x interval, slice size 
  unsigned nslice; // number of slices




public :
    //! -
  FastFinder(const BaseStarList &List);
  //! -
  const BaseStar *FindClosest(const Point &Where, const double MaxDist, 
			      bool (*SkipIt)(const BaseStar *) = NULL) const;


  const BaseStar *SecondClosest(const Point &Where, 
				const double MaxDist, 
				const BaseStar* &Closest,
  				bool (*SkipIt)(const BaseStar *)= NULL) const;



  void dump() const;  

  class Iterator {
    public : // could be made private, but what for??
    const FastFinder *finder;
    unsigned startSlice;
    unsigned endSlice;
    unsigned currentSlice;
    double yStart, yEnd;
    pstar current, pend, null_value;
    void check() const;

  public:
  Iterator(const FastFinder *f) : finder(f), null_value(f->stars.end()) {}
    void operator++() ;
    stars_element operator*() const;
  };

#ifndef SWIG
  Iterator begin_scan(const Point &Where, const double &MaxDist) const;
#endif

  void yslice(const int iSlice, const double YStart, const double YEnd, pstar &Start, pstar &End) const;
  pstar locate_y_start(pstar Begin, pstar End, const double &YVal) const;
  pstar locate_y_end(pstar Begin, pstar End, const double &YVal) const;

};

}}}
#endif /* FASTFINDER__H */
