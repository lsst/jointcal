#include <algorithm>

#include "lsst/meas/simastrom/BaseStar.h"
#include "lsst/meas/simastrom/FastFinder.h"

namespace lsst {
namespace meas {
namespace simastrom {



FastFinder::FastFinder(const BaseStarList &List, const unsigned NXslice) : baselist(List), count(List.size()), stars(count), nslice(NXslice), index(nslice+1)
{
  if (count==0) return;

  // fill "stars"
  unsigned j=0;
  for (BaseStarCIterator ci = List.begin(); ci != List.end(); ++ci)
    {
      stars[j] = ci->get();
      ++j;
    }

  sort(stars.begin(), stars.end(),
       [](const stars_element &E1, const stars_element &E2) 
       { return (E1->x < E2->x);} );

  xmin = stars[0]->x;
  xmax = stars[count-1]->x;
  nslice = std::min(nslice, count);
  if (xmin == xmax) nslice = 1;

  // the x size of each slice:
  xstep = (xmax-xmin)/nslice;
  

  // fill the index array with the first star beyond the slice limit.
  index[0] = 0; // first
  unsigned istar=0;
  for (unsigned islice=1; islice<nslice; ++islice)
    {
      double xend = xmin+(islice)*xstep;
      while (istar < count && stars[istar]->x < xend) ++istar;
      index[islice] = istar;
    }
  index[nslice] = count; // last
  for (unsigned islice=0; islice<nslice; ++islice)
    {
      sort(stars.begin()+index[islice], stars.begin()+index[islice+1], 
       [](const stars_element &E1, const stars_element &E2) 
       { return (E1->y < E2->y);} );// sort each slice in y.
    }
  //dump();
}


void FastFinder::dump() const
{
  for (unsigned i=0; i<count; ++i)
    {
      stars[i]->dump();
    }
}

const BaseStar *FastFinder::FindClosest(const Point &Where, 
					const double MaxDist, 
					bool (*SkipIt)(const BaseStar *)) const
{
  if (count == 0) return NULL;
  FastFinder::Iterator it = begin_scan(Where, MaxDist);
  if (*it == NULL) return NULL;
  const BaseStar *pbest = NULL;
  double minDist2 = MaxDist*MaxDist;
  for (      ; *it != NULL ; ++it)
    {
      const BaseStar *p = *it;
      if (SkipIt && SkipIt(p)) continue;
      double dist2 = Where.Dist2(*p);
      if (dist2 < minDist2) { pbest = p; minDist2 = dist2; } 
    }
  //  cout << "Distance " << minDist2 << " " << dist  << endl;
  return pbest;
}

const BaseStar *FastFinder::SecondClosest(const Point &Where, 
					  const double MaxDist, 
					  const BaseStar* &Closest,
					  bool (*SkipIt)(const BaseStar *)) const
{
  Closest=NULL;
  if (count == 0) return NULL;
  FastFinder::Iterator it = begin_scan(Where, MaxDist);
  if (*it == NULL) return NULL;
  const BaseStar *pbest1 = NULL; // closest
  const BaseStar *pbest2 = NULL; // second closest
  double minDist1_2 = MaxDist*MaxDist;
  double minDist2_2 = MaxDist*MaxDist;
  for (      ; *it != NULL ; ++it)
    {
      const BaseStar *p = *it;
      if (SkipIt && SkipIt(p)) continue;
      double dist2 = Where.Dist2(*p);
      if (dist2 < minDist1_2) 
	{ 
	  pbest2= pbest1;
	  minDist2_2 = minDist1_2;
	  pbest1 = p; 
	  minDist1_2 = dist2;
	}
      else if (dist2< minDist2_2)
	{
	  pbest2 = p;
	  minDist2_2 = dist2;
	}
    }
  Closest = pbest1;
  //  cout << "Distance " << minDist2 << " " << dist  << endl;
  return pbest2;
}

     

/* It is by no means clear the the 2 following routines are actually needed.
   It is nor clear to me (P.A) why they are different... but they really are.
*/
/* Locate the last position (in the sorted array) between Begin and
   End that lies before YVal.*/
FastFinder::pstar FastFinder::locate_y_start(pstar Begin, pstar End, const double &YVal) const
{
  if (Begin==stars.end() || Begin == End) return stars.end();
  int span = End - Begin -1;
  while (span > 1)
    {
      int half_span = span/2;
      pstar middle = Begin + half_span;
      if ((*middle)->y < YVal)
	{
	  Begin += half_span;
	  span -= half_span;
	}
      else
	{
	  span -= (span-half_span);
	}
    }
  return Begin;
}

/* Locate the first position (in the sorted array) between Begin and
   End that lies beyond YVal.*/
FastFinder::pstar FastFinder::locate_y_end(pstar Begin, pstar End, const double &YVal) const
{
  if (Begin==stars.end()) return stars.end();
  int span = End - Begin -1;
  while (span > 1)
    {
      int half_span = span/2;
      pstar middle = End - half_span;
      if ((*middle)->y > YVal)
	{
	  End -= half_span;
	  span -= half_span;
	}
      else
	{
	  span -= (span-half_span);
	}
    }
  return End-1;
}


void FastFinder::find_range_in_slice(const int iSlice, 
				     const double YStart, const double YEnd, 
				     pstar &Start, pstar &End) const
{
  Start = locate_y_start(stars.begin()+index[iSlice], stars.begin()+index[iSlice+1],  YStart);
  End   = locate_y_end(Start,stars.begin()+index[iSlice+1], YEnd);
}


FastFinder::Iterator  FastFinder::begin_scan(const Point &Where, const double &MaxDist) const
{
  return FastFinder::Iterator(*this,Where, MaxDist);
}

using Iterator = FastFinder::Iterator;


Iterator::Iterator(const FastFinder &F, const Point &Where, 
		   const double &MaxDist) 
  : finder(F), null_value(F.stars.end()) 
{
  current = pend = null_value;// does not iterate
  int startSlice = 0;
  if (finder.xstep != 0) // means we have several slices
    {
      startSlice = std::max(0,int((Where.x-MaxDist-finder.xmin)/finder.xstep));
      /* obviously, endSlice (and starSlice) can be negative. 
	 This is why slice indices are "int" rather than "unsigned". */
      endSlice = std::min(int(finder.nslice),int((Where.x+MaxDist-finder.xmin)/finder.xstep)+1);
    }
  else 
    {
      startSlice = 0;
      endSlice = 1;
    }
  // beyond limits: 
  if (startSlice >= int(finder.nslice) || endSlice < 0) return;
  // we are inside in x, so, we setup the y range:
  yStart = Where.y - MaxDist;
  yEnd   = Where.y + MaxDist;
  /* rather than initializing here, we step back one 
     slice and let "++" do its job */
  currentSlice = startSlice -1 ; // again, this requires "int" slices
  ++(*this);
}

FastFinder::stars_element  Iterator::operator*() const
{
  if (current!=null_value) return *current; 
  return NULL;
}

void Iterator::operator++()
{
  if (current != pend) {current++;}
  else
  do 
    {
      currentSlice++;
      if (currentSlice >= endSlice) {current = null_value; return;}
      finder.find_range_in_slice(currentSlice, yStart, yEnd, current, pend);
    } while(current == null_value);
  check();
}


void FastFinder::Iterator::check() const
{
  if (current != null_value && (current < finder.stars.begin() || current >= finder.stars.begin()+finder.count))
    {
      std::cout << "ERROR in FastFinder " << *current << " " << *(finder.stars.begin()) << ' ' << *(finder.stars.begin()+finder.count) << std::endl;
    }
}

}}} // end of namespaces
