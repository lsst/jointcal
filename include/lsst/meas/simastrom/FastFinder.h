#ifndef FASTFINDER__H
#define FASTFINDER__H

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
  const BaseStar **stars;        // array of pointers to objects in the StarList 
  double xmin,xmax, xstep; // x interval, slice size 
  int nslice; // number of slices
  int *index;              // index in "stars" of first object of each slice.
  int count;               // total number of objects (size of stars).
  const BaseStarList* baselist; // pointer to the initial list of stars

public :
    //! -
  FastFinder(const BaseStarList &List);
  ~FastFinder();
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
    int startSlice;
    int endSlice;
    int currentSlice;
    double yStart, yEnd;
    const BaseStar **current, **pend;
    const FastFinder *finder;
    void check() const;

  public:
    void operator++() ;
    const BaseStar* operator*() const;
  };

#ifndef SWIG
  Iterator begin_scan(const Point &Where, const double &MaxDist) const;
#endif

  void yslice(const int iSlice, const double YStart, const double YEnd, const BaseStar **&Start, const BaseStar **&End) const;
  const BaseStar** locate_y_start(const BaseStar **Begin, const BaseStar **End, const double &YVal) const;
  const BaseStar** locate_y_end(const BaseStar **Begin, const BaseStar **End, const double &YVal) const;

};

}}}
#endif /* FASTFINDER__H */
