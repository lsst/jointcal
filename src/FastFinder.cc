#include <algorithm>

#include "lsst/meas/simastrom/BaseStar.h"
#include "lsst/meas/simastrom/FastFinder.h"

namespace lsst {
namespace meas {
namespace simastrom {


using namespace std;

// could become an argument of the constructor..
#define NSLICE 100U


// used to sort the array of pointers.
static bool CompareX(const BaseStar* p1, const BaseStar* p2)
{
  return (p1->x < p2->x);
}


static bool CompareY(const BaseStar* p1, const BaseStar* p2)
{
  return (p1->y < p2->y);
}


FastFinder::FastFinder(const BaseStarList &List) : baselist(List), count(List.size()), stars(count), index(nslice+1)
{
  // fill "stars"
  int j=0;
  for (BaseStarCIterator ci = List.begin(); ci != List.end(); ++ci)
    {
      stars[j] = ci->get();
      j++;
    }

  sort(stars.begin(), stars.end(),[](const stars_element &E1, const stars_element &E2) -> bool { return (E1->x < E2->x);} );

  xmin = stars[0]->x;
  xmax = stars[count-1]->x;
  nslice = min(NSLICE,count);
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
      sort(stars.begin()+index[islice], stars.begin()+index[islice+1], CompareY);  // sort each slice in y.
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


void FastFinder::yslice(const int iSlice, const double YStart, const double YEnd, pstar &Start, pstar &End) const
{
  Start = locate_y_start(stars.begin()+index[iSlice], stars.begin()+index[iSlice+1],  YStart);
  End   = locate_y_end(Start,stars.begin()+index[iSlice+1], YEnd);
}

FastFinder::Iterator  FastFinder::begin_scan(const Point &Where, const double &MaxDist) const
{
  FastFinder::Iterator iterator(this);
  if (xstep != 0)
    {
      iterator.startSlice = max(0,int((Where.x-MaxDist-xmin)/xstep));
      int endslice = min(int(nslice),int((Where.x+MaxDist-xmin)/xstep)+1);
      iterator.endSlice = (endslice > 0) ? endslice : 0;
    }
  else 
    {
      iterator.startSlice = 0;
      iterator.endSlice = 1;
    }
  iterator.current = stars.end();
  if (iterator.startSlice >= nslice || iterator.endSlice < 0) return iterator;
  iterator.current = stars.end();
  iterator.yStart = Where.y - MaxDist;
  iterator.yEnd   = Where.y + MaxDist;
  iterator.current = iterator.pend = stars.end();
  iterator.currentSlice = iterator.startSlice - 1;
  ++iterator;
  return iterator;
}

FastFinder::stars_element  FastFinder::Iterator::operator*() const
{
  if (current==null_value) return *current; else return NULL;
}

void FastFinder::Iterator::operator++()
{
  if (current != pend) {current++;}
  else
  do 
    {
      currentSlice++;
      if (currentSlice >= endSlice) {current = null_value; return;}
      finder->yslice(currentSlice, yStart, yEnd, current, pend);
    } while(current == null_value);
  check();
}


void FastFinder::Iterator::check() const
{
  if (current != null_value && (current < finder->stars.begin() || current >= finder->stars.begin()+finder->count))
    {
      std::cout << " alerte !! " << *current << " " << *(finder->stars.begin()) << ' ' << *(finder->stars.begin()+finder->count) << std::endl;
      exit(1);
    }
}

}}} // end of namespaces
