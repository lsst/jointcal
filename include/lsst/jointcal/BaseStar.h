// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_BASE_STAR_H
#define LSST_JOINTCAL_BASE_STAR_H


#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>

#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {


#define MEMPIX2DISK 1


#define DECALAGE_IJ_XY 0.
#define DECALAGE_XY_IJ 0.
/*! \file */

// tell other pieces of code that BaseStar now derives from FatPoint (instead of Point)
#define BASESTAR_HAS_POSITION_ERRORS

//! The base class for handling stars. Used by all matching routines.
class BaseStar : public FatPoint
{

  public :
double flux;


  public:
 BaseStar(){x=0;y=0;flux=0;};
  //! constructor
 BaseStar(double xx, double yy, double ff) : FatPoint(xx,yy), flux(ff)
  {};
 BaseStar(const Point &a_point, double a_flux) : FatPoint(a_point), flux(a_flux)
  {};

  //! access stuff.
  double X() const { return x;}
  //!
  double Y() const { return y;}

  //! allows std::cout << aBaseStar;
  friend std::ostream& operator << (std::ostream &stream, const BaseStar &s)
  { s.dump(stream); return stream;}

  virtual std::string __str__() const {std::stringstream s; dump(s); return s.str();}

  virtual void dump(std::ostream & stream = std::cout) const { stream << "x: "<< x << " y: " << y << " flux: " << flux;}

  BaseStar& operator=(const Point &P) {this->x = P.x; this->y = P.y; return (*this);};

  static const char *TypeName() { return "BaseStar";}

  virtual ~BaseStar(){};
};


//! enables to sort easily a starList (of anything that derives from BaseStar)
bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2);

int DecodeFormat(const char *FormatLine, const char *StarName);


typedef StarList<BaseStar> BaseStarList;

typedef BaseStarList::const_iterator BaseStarCIterator;
 typedef BaseStarList::iterator BaseStarIterator;


}}

#endif // LSST_JOINTCAL_BASE_STAR_H
