// -*- C++ -*-
#ifndef REFSTAR__H
#define REFSTAR__H

#include <vector>
#include <fstream>

#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {


/*! \file */

//! Objects used as position anchors, typically USNO stars. Coordinate system defined by user. The Common Tangent Plane seems a good idea.
class RefStar : public BaseStar
{
private :
  unsigned int index;
  FittedStarRef fittedStar;
  Point raDec;
  std::vector<double> refFlux;

  public :

  RefStar();

  //!
  RefStar(const BaseStar &, const Point &RaDec);

  //!
  void SetFittedStar(FittedStar &F);

  //!
  double Ra() const { return raDec.x;}

  //!
  double Dec() const {return raDec.y;}

  //! reference flux
  double Flux(int band) const;

  //! assign the reference fluxes
  void   AssignRefFluxes(std::vector<double> const& refflux);

  //! star index
  unsigned int &  Index() { return index; }
  unsigned int    Index() const { return index; }

};



class RefStarList : public StarList<RefStar> { };

typedef RefStarList::const_iterator RefStarCIterator;
typedef RefStarList::iterator RefStarIterator;
typedef CountedRef<RefStar> RefStarRef;



/********* RefStarTuple ****************/

class RefStarTuple {
  private :
  std::ofstream stream;

  public :
   RefStarTuple( const std::string &FileName);

  void AddEntry(const RefStar &R, const FittedStar &F);

  ~RefStarTuple() { stream.close();}
};

}} // end of namespaces

#endif /* REFSTAR__H */
