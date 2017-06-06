// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_REF_STAR_H
#define LSST_JOINTCAL_REF_STAR_H

#include <vector>
#include <fstream>

#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {

//! Objects used as position anchors, typically USNO stars. Coordinate system defined by user. The Common
//! Tangent Plane seems a good idea.
class RefStar : public BaseStar {
private:
    unsigned int _index;
    std::vector<double> _refFlux;

public:
    //!
    RefStar(const BaseStar& baseStar);

    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " refFlux: [";
        for (auto x : _refFlux) {
            stream << x << ", ";
        }
        stream << "] index: " << _index;
    }

    //! reference flux
    double getFlux(int filter) const;

    //! assign the reference fluxes
    void assignRefFluxes(std::vector<double> const& refFlux);

    //! star index
    unsigned int& getIndex() { return _index; }
    unsigned int getIndex() const { return _index; }
};

/****** RefStarList ***********/

class Frame;

// typedef StarList<RefStar> RefStarList;
class RefStarList : public StarList<RefStar> {};

typedef RefStarList::const_iterator RefStarCIterator;
typedef RefStarList::iterator RefStarIterator;

BaseStarList& Ref2Base(RefStarList& This);
BaseStarList* Ref2Base(RefStarList* This);
const BaseStarList& Ref2Base(const RefStarList& This);
const BaseStarList* Ref2Base(const RefStarList* This);
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_REF_STAR_H
