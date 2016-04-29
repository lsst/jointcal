#ifndef CCDIMAGELIST__H
#define CCDIMAGELIST__H

#include <list>

#include "lsst/jointcal/CcdImage.h"


namespace lsst {
namespace jointcal {

//! a  list of CcdImage. Usually produced by Associations
//class CcdImageList : public std::list<CountedRef<CcdImage> >
class CcdImageList : public std::list<boost::shared_ptr<CcdImage> >
{
  public:

  //!
  //std::list<std::string> DateObs() const;

  //!
  //std::list<std::string> Bands() const;

  //!
  //double       MeanAirmass() const;

  //!
  template<class Accept> CcdImageList SubList(const Accept &OP) const
    {
      CcdImageList out;
      for (const_iterator i = begin(); i != end() ; ++i)
	if (OP(**i)) out.push_back(*i);
      return out;
    }

  // find the matching image. Chip==-1 means any chip
//  double AirMass(const int Shoot, const int Chip = -1) const;
};


typedef CcdImageList::iterator CcdImageIterator;
typedef CcdImageList::const_iterator CcdImageCIterator;

}} // end of namespaces

#endif /* CCDIMAGELIST__H */
