// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_LIST_MATCH_H
#define LSST_JOINTCAL_LIST_MATCH_H

#include <string>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarMatch.h"

namespace lsst {
namespace jointcal {

class Gtransfo;
class GtransfoLin;

//! Parameters to be provided to combinatorial searches
struct MatchConditions {
    int NStarsL1, NStarsL2;
    int MaxTrialCount;
    double NSigmas;
    double MaxShiftX, MaxShiftY;
    double SizeRatio, DeltaSizeRatio, MinMatchRatio;
    int PrintLevel;
    int Algorithm;

    MatchConditions(/* const std::string &DatacardsName = ""*/);

    double MinSizeRatio() const { return SizeRatio - DeltaSizeRatio; }
    double MaxSizeRatio() const { return SizeRatio + DeltaSizeRatio; }
};

/*! \file
    \brief Combinatorial searches for linear transformations to go from
           list1 to list2.

    The following routines search a geometrical transformation that make
two lists of stars to match geometrically as well as possible. They are used
either to match two images of the same sky area, or an image with a catalogue.
They assume that fluxes assigned to stars are actual fluxes, i.e. the brighter
the star, the higher the flux. They however only rely on flux ordering,
not values.
 */

//! searches a geometrical transformation that goes from List1 to List2.
/*!  The found transformation is a field of the returned object, as well as the star pairs
(the matches) that were constructed.  (see StarMatchList class definition for more details).
The various cuts are contained in Conditions (see listmatch.h) for its contents.
This routine searches a transformation that involves a shift and a rotation. */

std::unique_ptr<StarMatchList> MatchSearchRotShift(BaseStarList &L1, BaseStarList &L2,
                                                   const MatchConditions &Conditions);

//! same as above but searches also a flipped solution.

std::unique_ptr<StarMatchList> MatchSearchRotShiftFlip(BaseStarList &L1, BaseStarList &L2,
                                                       const MatchConditions &Conditions);

//! assembles star matches.
/*! It picks stars in L1, transforms them through Guess, and collects
closest star in L2, and builds a match if closer than MaxDist). */

std::unique_ptr<StarMatchList> ListMatchCollect(const BaseStarList &L1, const BaseStarList &L2,
                                                const Gtransfo *Guess, const double MaxDist);

//! same as before except that the transfo is the identity

std::unique_ptr<StarMatchList> ListMatchCollect(const BaseStarList &L1, const BaseStarList &L2,
                                                const double MaxDist);

//! searches for a 2 dimensional shift using a very crude histogram method.

std::unique_ptr<GtransfoLin> ListMatchupShift(const BaseStarList &L1, const BaseStarList &L2,
                                              const Gtransfo &Tin, double MaxShift, double BinSize = 0);

std::unique_ptr<Gtransfo> ListMatchCombinatorial(const BaseStarList &List1, const BaseStarList &List2,
                                                 const MatchConditions &Conditions = MatchConditions());
std::unique_ptr<Gtransfo> ListMatchRefine(const BaseStarList &List1, const BaseStarList &List2,
                                          std::unique_ptr<Gtransfo> transfo, const int maxOrder = 3);

#ifdef DO_WE_NEED_THAT
inline Gtransfo *ListMatch(const BaseStarList &List1, const BaseStarList &List2, const int maxOrder = 3) {
    Gtransfo *transfo = ListMatchCombinatorial(List1, List2);
    transfo = ListMatchRefine(List1, List2, transfo, maxOrder);
    return transfo;
}
#endif /*  DO_WE_NEED_THAT */
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_LIST_MATCH_H
