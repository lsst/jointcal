// This may look like C code, but it is really -*- C++ -*-
#ifndef SIMASTROM_POINT__H
#define SIMASTROM_POINT__H
#include <iostream>
#include <cmath>

namespace lsst {
namespace meas {
namespace simastrom {


/*! \file */

//! A point in a plane.
class Point {
public:
  virtual ~Point() {}

  //! coordinate
  double x,y;

  //! - contructor
  Point() : x(0), y(0) {};

  //! - contructor
  Point(const double &xx, const double &yy) : x(xx), y(yy) {};


  //! -
   double Distance(const Point& Other) const { return sqrt((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! distance squared to Other
   double Dist2(const Point& Other) const { return ((x-Other.x)*(x-Other.x) + (y-Other.y)*(y-Other.y));};

  //! Sum
  Point operator + (const Point &Right) const { return Point(x+Right.x, y+Right.y);}


  //! Difference
  Point operator - (const Point &Right) const { return Point(x-Right.x, y-Right.y);}

  //! utility
  virtual void dump(std::ostream& s = std::cout) const { s <<" x " << x << " y " << y;}  

#ifndef SWIG
  //! -
  friend std::ostream& operator << (std::ostream& stream, const Point &point) 
  { point.dump(stream); return stream;}
#endif
};

#ifndef SWIG
//#ifdef __CINT__ /* rootcint  does not recognise friend inline functions */
std::ostream& operator << (std::ostream& stream, const Point &point);
//#endif
#endif

}}}

#endif /* POINT__H */


