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


template<class Star>void StarList<Star>::FluxSort()
{
  typedef StarList<Star>::Element E;
  this->sort([] (const E &E1, const E &E2) {return (E1->flux > E2->flux);});
}

template<class Star>void StarList<Star>::CutTail(const int NKeep)
{
  int count = 0;
  auto si = this->begin();
  for (  ; si != this->end() && count < NKeep; ++count, ++si);
  while (si != this->end() ) {si = this->erase(si);}
}

template<class Star>void StarList<Star>::ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const
{
  for (auto const &star: *this)
    {
      if (aFrame.InFrame(*star))
	{
	  Star *copy = new Star(*star);
	  Out.push_back(copy);
	}
    }
}

template<class Star>void StarList<Star>::CopyTo(StarList<Star> &Copy) const
{
  Copy.ClearList();
  //  Copy.GlobVal() = this->GlobVal();
  for (auto const &si: *this)
    Copy.push_back(new Star(*si));
}

// Explicit instantiations
template class StarList<BaseStar>;
template class StarList<FittedStar>;
template class StarList<MeasuredStar>;

}} // end of namespaces

#endif /* STARLIST__CC */
