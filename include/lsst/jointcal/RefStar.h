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
    // on-sky flux, in Maggies, per filter
    std::vector<double> _refFluxList;
    std::vector<double> _refFluxErrList;

public:
    // RefStar(double xx, double yy, std::vector<double>& refFluxList, std::vector<double>& refFluxErrList)
    //         : FatPoint(xx, yy), _refFluxList(refFluxList), _refFluxErrList(refFluxErrList) {}

    RefStar(double xx, double yy, double defaultFlux, std::vector<double>& refFluxList,
            std::vector<double>& refFluxErrList)
            : BaseStar(xx, yy, defaultFlux), _refFluxList(refFluxList), _refFluxErrList(refFluxErrList) {}
    //!
    RefStar(const BaseStar& baseStar) : BaseStar(baseStar), _index(0) {}

    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " refFlux: [";
        for (auto x : _refFluxList) {
            stream << x << ", ";
        }
        stream << "] index: " << _index;
    }

    /// reference flux in a given filter
    double getFlux(size_t filter) const { return _refFluxList[filter]; }
    /// reference fluxErr in a given filter
    double getFluxErr(size_t filter) const { return _refFluxErrList[filter]; }

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
