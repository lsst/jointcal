#ifndef STARLIST__CC
#define STARLIST__CC

#include "lsst/pex/exceptions.h"

#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

template <class Star>
void StarList<Star>::fluxSort() {
    typedef StarList<Star>::Element E;
    this->sort([](const E &e1, const E &e2) { return (e1->getFlux() > e2->getFlux()); });
}

template <class Star>
void StarList<Star>::cutTail(const int nKeep) {
    int count = 0;
    auto si = this->begin();
    for (; si != this->end() && count < nKeep; ++count, ++si)
        ;
    while (si != this->end()) {
        si = this->erase(si);
    }
}

template <class Star>
void StarList<Star>::extractInFrame(StarList<Star> &out, const Frame &frame) const {
    for (auto const &star : *this) {
        if (frame.inFrame(*star)) {
            out.push_back(std::make_shared<Star>(*star));
        }
    }
}

template <class Star>
void StarList<Star>::copyTo(StarList<Star> &copy) const {
    copy.clearList();
    for (auto const &si : *this) copy.push_back(std::make_shared<Star>(*si));
}

// Explicit instantiations
template class StarList<BaseStar>;
template class StarList<FittedStar>;
template class StarList<MeasuredStar>;
}  // namespace jointcal
}  // namespace lsst

#endif /* STARLIST__CC */
