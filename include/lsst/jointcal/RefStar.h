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
  std::vector<double> refFlux;

  public :
  //!
  RefStar(const BaseStar &baseStar);

  void dump(std::ostream & stream = std::cout) const
    {
        BaseStar::dump(stream);
        stream << " refFlux: [";
        for (auto x: refFlux) {
            stream << x << ", ";
        }
        stream << "] index: " << index;
    }

  //! reference flux
  double Flux(int filter) const;

  //! assign the reference fluxes
  void   AssignRefFluxes(std::vector<double> const& refflux);

  //! star index
  unsigned int &  Index() { return index; }
  unsigned int    Index() const { return index; }

};

/****** RefStarList ***********/


class Frame;

  //typedef StarList<RefStar> RefStarList;
class RefStarList : public StarList<RefStar> { };


typedef RefStarList::const_iterator RefStarCIterator;
typedef RefStarList::iterator RefStarIterator;


BaseStarList& Ref2Base(RefStarList &This);
BaseStarList* Ref2Base(RefStarList *This);
const BaseStarList& Ref2Base(const RefStarList &This);
const BaseStarList* Ref2Base(const RefStarList *This);

}} // end of namespaces

#endif /* REFSTAR__H */
