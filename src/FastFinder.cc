#include <algorithm>

#include "lsst/meas/simastrom/BaseStar.h"
#include "lsst/meas/simastrom/FastFinder.h"

namespace lsst {
namespace meas {
namespace simastrom {


using namespace std;

// could become an argument of the constructor..
#define NSLICE 100


// used to sort the array of pointers.
static bool CompareX(const BaseStar* p1, const BaseStar* p2)
{
  return (p1->x < p2->x);
}


static bool CompareY(const BaseStar* p1, const BaseStar* p2)
{
  return (p1->y < p2->y);
}


FastFinder::FastFinder(const BaseStarList &List)
{
  baselist = &List;
  count = List.size();
  stars = NULL;
  index = NULL;
  if (count == 0) return;

  stars = new const BaseStar*[count];

  // fill "stars"
  int j=0;
  for (BaseStarCIterator ci = List.begin(); ci != List.end(); ++ci)
    {
      stars[j] = *ci;
      j++;
    }
  sort(stars,stars+count,CompareX);

  xmin = stars[0]->x;
  xmax = stars[count-1]->x;
  nslice = min(NSLICE,count);
  if (xmin == xmax) nslice = 1;

  // index enables fast access to a given x slice.
  index = new int[nslice+1]; // "+1" because index contains the limits of NSLICE slices.

  // the x size of each slice:
  xstep = (xmax-xmin)/nslice;
  

  // fill the index array with the first star beyond the slice limit.
  index[0] = 0; // first
  int istar=0;
  for (int islice=1; islice<nslice; ++islice)
    {
      double xend = xmin+(islice)*xstep;
      while (istar < count && stars[istar]->x < xend) ++istar;
      index[islice] = istar;
    }
  index[nslice] = count; // last
  for (int islice=0; islice<nslice; ++islice)
    {
      sort(stars+index[islice], stars+index[islice+1], CompareY);  // sort each slice in y.
    }
  //dump();
}

FastFinder::~FastFinder()
{ 
  if (stars) delete [] stars; 
  if (index) delete [] index;}


void FastFinder::dump() const
{
  for (int i=0; i<count; ++i)
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

const BaseStar **FastFinder::locate_y_start(const BaseStar **Begin, const BaseStar **End, const double &YVal) const
{
  if (!Begin || Begin == End) return NULL;
  int span = End - Begin -1;
  while (span > 1)
    {
      int half_span = span/2;
      const BaseStar **middle = Begin + half_span;
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


const BaseStar **FastFinder::locate_y_end(const BaseStar **Begin, const BaseStar **End, const double &YVal) const
{
  if (!Begin) return NULL;
  int span = End - Begin -1;
  while (span > 1)
    {
      int half_span = span/2;
      const BaseStar **middle = End - half_span;
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


void FastFinder::yslice(const int iSlice, const double YStart, const double YEnd, const BaseStar **&Start, const BaseStar **&End) const
{
  Start = locate_y_start(stars+index[iSlice], stars+index[iSlice+1],  YStart);
  End   = locate_y_end(Start,stars+index[iSlice+1], YEnd);
}

FastFinder::Iterator  FastFinder::begin_scan(const Point &Where, const double &MaxDist) const
{
  FastFinder::Iterator iterator;
  iterator.finder = this;
  if (xstep != 0)
    {
      iterator.startSlice = max(0,int((Where.x-MaxDist-xmin)/xstep));
      iterator.endSlice = min(nslice,int((Where.x+MaxDist-xmin)/xstep)+1);
    }
  else 
    {
      iterator.startSlice = 0;
      iterator.endSlice = 1;
    }
  iterator.current = NULL;
  if (iterator.startSlice >= nslice || iterator.endSlice < 0) return iterator;
  iterator.current = NULL;
  iterator.yStart = Where.y - MaxDist;
  iterator.yEnd   = Where.y + MaxDist;
  iterator.current = iterator.pend = NULL;
  iterator.currentSlice = iterator.startSlice - 1;
  ++iterator;
  return iterator;
}

const BaseStar*  FastFinder::Iterator::operator*() const
{
  if (current) return *current; else return NULL;
}

void FastFinder::Iterator::operator++()
{
  if (current != pend) {current++;}
  else
  do 
    {
      currentSlice++;
      if (currentSlice >= endSlice) {current = NULL; return;}
      finder->yslice(currentSlice, yStart, yEnd, current, pend);
    } while(current == NULL);
  check();
}


void FastFinder::Iterator::check() const
{
  if (current && (current < finder->stars || current >= finder->stars+finder->count))
    {
      std::cout << " alerte !! " << current << " " << finder->stars << ' ' << finder->stars+finder->count << std::endl;
    }
}

}}} // end of namespaces
