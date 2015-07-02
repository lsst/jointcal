
// This may look like C code, but it is really -*- C++ -*-


#ifndef BASESTAR__H
#define BASESTAR__H


#include <iostream>
#include <cstdio>
#include <string>

#include "lsst/meas/simastrom/FatPoint.h"
#include "lsst/meas/simastrom/CountedRef.h"

namespace lsst {
namespace meas {
namespace simastrom {


#define MEMPIX2DISK 1


#define DECALAGE_IJ_XY 0.
#define DECALAGE_XY_IJ 0.
/*! \file */

// tell other pieces of code that BaseStar now derives from FatPoint (instead of Point)
#define BASESTAR_HAS_POSITION_ERRORS

//! The base class for handling stars. Used by all matching routines.
class BaseStar : public FatPoint, public RefCount
{

  public : // si quelqu'un connait un moyen efficace d'eviter ca...
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

#ifdef DO_WE_NEED_IT
  static BaseStar* read(fastifstream & rd, const char *format);

  void  read_it(fastifstream & rd, const char *format);
#endif /* DO_WE_NEED_IT */

  virtual void write(std::ostream &s = std::cout)const ;
  virtual void writen(std::ostream &s = std::cout)const ;

#ifndef SWIG
  //! allows std::cout << aBaseStar;
  friend std::ostream& operator << (std::ostream &stream, const BaseStar &s)
  { s.dump(stream); return stream;}
#endif

  virtual void dump(std::ostream & stream = std::cout) const { stream << "x "<< x << " y " << y << " flux " << flux << std::endl;}
  virtual void dumpn(std::ostream & stream = std::cout) const { stream << "x "<< x << " y " << y << " flux " << flux << " ";}

#ifndef SWIG
  BaseStar& operator=(const Point &P) {this->x = P.x; this->y = P.y; return (*this);};
#endif

  static const char *TypeName() { return "BaseStar";}

  virtual ~BaseStar(){};

  virtual std::string WriteHeader_(std::ostream & stream = std::cout, const char*i = NULL) const ;

  virtual void WriteHeader(std::ostream & stream = std::cout) const;
#ifdef USE_ROOT

#endif
  //  ClassDef(BaseStar,1) // no ";" ....
};

//! Number of values read for this format
unsigned NValsBaseStar(const char *Format);


//! enables to sort easily a starList (of anything that derives from BaseStar) 
bool DecreasingFlux(const BaseStar *S1, const BaseStar *S2);



#ifdef TO_BE_FIXED
int DecodeFormat(const char *FormatLine, const char *StarName);


/* what concerns the BaseStarList's : */
#include "starlist.h"

typedef StarList<BaseStar> BaseStarList;



typedef BaseStarList::const_iterator BaseStarCIterator;
typedef BaseStarList::iterator BaseStarIterator;
typedef CountedRef<BaseStar> BaseStarRef;

#ifdef __GCCXML__
template class CountedRef<BaseStar>;
template class std::list<BaseStar*>;
// template class std::list<CountedRef<BaseStar> >;
#endif
#endif /* TO_BE_FIXED */

}}}

#endif /* BASESTAR__H */  
