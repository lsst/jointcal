// This may look like C code, but it is really -*- C++ -*-
#ifndef FRAME__H
#define FRAME__H

#include <iostream>
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {


typedef enum {WholeSizeFrame, ClippedSizeFrame} WhichFrame;


typedef enum {LargeFrame, SmallFrame} WhichTransformed;

//! rectangle with sides parallel to axes.
/*! when Frame's are used to define subparts of images, xMin and xMax refer
  to the first and last pixels in the subimage */

class Frame {
public:
  //! coordinate of boundary.
  double xMin,xMax,yMin,yMax;
  
  //! Default constructor
  Frame();
    
  //! this one is dangerous: you may swap the 2 middle arguments. 
  //! Prefer next one
  Frame(const double &xMin, const double &yMin, 
	const double &xMax, const double &yMax);
  
  //! typical use: Frame(Point(xmin,ymin),Point(xmax,ymax))
  Frame(const Point &LowerLeft, const Point &UpperRight);

#ifdef TO_BE_FIXED
  //! 2 kinds of bounds in headers, the chip size and some 
  //! that may be added by hand. See WriteInHeader().
  Frame(const FitsHeader &header, WhichFrame which=ClippedSizeFrame);
#endif

  //! size along x axis
  double Width() const {return xMax-xMin;}

  //! size along y axis
  double Height() const {return yMax-yMin;}
  
  //! Center  of the frame
  Point Center() const {return Point((xMax+xMin)*0.5,(yMax+yMin)*0.5);}
  
  //! intersection of Frame's.
  Frame operator*(const Frame &Right) const;  /* intersection : a = b n c */
  
  //! intersection of Frame's
  Frame& operator*=( const Frame &Right);     /* intersection : a = a n b */
  
  //! union of Frames
  Frame operator+(const Frame &Right) const;  /* union : a = b u c */

  //! union of Frames
  Frame& operator+=( const Frame &Right);     /* intersection : a = a u b */
  
  //! shrinks the frame (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const double MarginSize);
  
  //! shrinks the frame (if MarginSize>0), enlarges it (if MarginSize<0).
  void CutMargin(const double MarginX, const double MarginY);
  
  //! necessary for comparisons (!= is defined from this one implicitely)
  bool operator ==(const Frame &Right) const;
  
  //! comparison
  bool operator !=(const Frame &Right) const {return !(*this == Right);}

  //! rescale it. The center does not move.
  Frame Rescale(const double Factor) const;
  
  // the area.
  double Area() const;
  
  //! inside?
  bool InFrame(const double &x, const double &y) const;
  
  //! same as above
  bool InFrame(const Point &pt) const 
  {return InFrame(pt.x,pt.y);}
  
  //! distance to closest boundary.
  double MinDistToEdges(const Point &P) const;
  
  void dump(std::ostream & stream = std::cout) const;
  
  //! allows \verbatim std::cout << frame; \endverbatim.
  friend std::ostream & operator<<(std::ostream &stream, const Frame &Right) 
          { Right.dump(stream); return stream;};
  
private:
  void order();
};

}} // end of namespaces

#endif /* FRAME__H */
