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
    // on-sky flux, in Maggies, per filter
    std::vector<double> _refFluxList;
    std::vector<double> _refFluxErrList;

public:
    RefStar(double xx, double yy, double defaultFlux, double defaultFluxErr, std::vector<double>& refFluxList,
            std::vector<double>& refFluxErrList)
            : BaseStar(xx, yy, defaultFlux, defaultFluxErr),
              _refFluxList(refFluxList),
              _refFluxErrList(refFluxErrList) {}

    /// No move or copy: each RefStar is unique, and should be accessed/managed via shared_ptr.
    RefStar(RefStar const&) = delete;
    RefStar(RefStar&&) = delete;
    RefStar& operator=(RefStar const&) = default;
    RefStar& operator=(RefStar&&) = delete;

    void dump(std::ostream& stream = std::cout) const {
        BaseStar::dump(stream);
        stream << " refFlux: [";
        for (auto x : _refFluxList) {
            stream << x << ", ";
        }
        stream << "]";
    }

    using BaseStar::getFlux;
    using BaseStar::getFluxErr;
    /// reference flux in a given filter
    double getFlux(size_t filter) const { return _refFluxList[filter]; }
    /// reference fluxErr in a given filter
    double getFluxErr(size_t filter) const { return _refFluxErrList[filter]; }
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
