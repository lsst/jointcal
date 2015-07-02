#ifndef STARLIST__CC
#define STARLIST__CC
#include <iostream>
#include <iomanip>
#include <istream>
#include <string>
#include <cstdio>


#include "lsst/meas/simastrom/StarList.h"
#include "lsst/meas/simastrom/Frame.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace meas {
namespace simastrom {

//#include "starlistexception.h"

template<class Star> StarList<Star>::StarList(const std::string &FileName) /* to be changed if switch to Star rather than pointers to Stars */ 
{
  this->read(FileName);
}


template<class Star> int StarList<Star>::read(std::istream & r)
{
  char c ;
  char buff[4096];
  char format[4096];
  ClearList();
  while( r >> c ) // to test eof
    {
      r.unget() ;
      if ( (c == '@') ) 
	{
	  r.getline(buff,4096); 
	  glob.ProcessLine(buff);
	  continue;
	}
      if ( (c == '#') ) // we jump over the line  (not always ...)
        {
	  r.getline(buff,4096);
	  /* hack something reading " format <StarType> <integer>" to drive the decoding (in Star::read) */
	  char *p = strstr(buff,"format");
	  if (p) /* this test is enough because the format is the last line of the header ... */
          {
	    strncpy(format,p + strlen("format"),4096);
          }
        }
      else // actual star data
	{
	  Star* s = dynamic_cast<Star*>(Star::read(r, format));
	  if (!r) // means I/O issue 
	    {
	      delete s;
	      return 0;
	    }
	  if (!s) return 0; // other kind of I/O issue. 
	  // read the end of line, in case there are some items left
	  r.getline(buff,4096);
	  push_back(s); // append just read star to the std::list.
	}
    }
  return 1;
}

template<class Star> int 
StarList<Star>::read(const std::string &FileName)
{
  return ascii_read(FileName);
}

template<class Star> int 
StarList<Star>::ascii_read(const std::string &FileName)
{
  std::istream rd(FileName.c_str());
  if (!rd)
    {
      std::cout << " Starlist : cannot open " << FileName << std::endl;
      throw LSST_EXCEPT(pexExcept::NotFoundError, "StarList :cannot open file="+FileName);
    }
  int count = read(rd);
  if (rd.fail())
    throw LSST_EXCEPT(pexExcept::IoError, "bad extraction in StarList reader, file="+FileName);
  return count;
}

template<class Star> int StarList<Star>::write(const std::string &FileName) const
{
  std::ofstream pr(FileName.c_str());
  if (!pr)
    {
      std::cerr << " StarList::write : could not open " << FileName << std::endl;
      return 0;
    }
  write(pr);
  pr.close();
  return 1;
}



template<class Star> int StarList<Star>::write(std::ostream & pr) const
{
  std::ios::fmtflags  old_flags =  pr.flags(); 
  pr  << resetiosflags(std::ios::scientific) ;
  // pr  << setiosflags(0) ;
  pr  << setiosflags(std::ios::fixed) ;
  int oldprec = pr.precision();
  pr<< std::setprecision(8);

  // write GlobalValues if any
  std::vector<std::string> globs = glob.OutputLines();
  for (unsigned k = 0; k < globs.size(); ++k)
    pr << '@' << globs[k] << std::endl;

  // cannot use front() to detect emptyness
  if (this->begin() == this->end()) // empty std::list, and faster than (size() == 0)
    {
      Star dummy;
      dummy.WriteHeader(pr);
    }
  else this->front()->WriteHeader(pr);
  for (StarCIterator it= this->begin(); it!=this->end() ; it++ )
    {    
      (*it)->write(pr);
    }
  pr.flags(old_flags);
  pr << std::setprecision(oldprec);
 return 1;
}






template <class Star>  void StarList<Star>::ExtractHead(StarList<Star> &Out, int NHead) const
{
int count = 0;
for (StarCIterator s= this->begin(); s!= this->end() && (count < NHead); ++s, count++)
  {
  Star *copy = new Star(*(*s));  /* to be changed if switch to Star rather than pointers to Stars */
  Out.push_back(copy);
  }
}

template<class Star> Star* StarList<Star>::FindClosest(double X, double Y) const
{
double min_dist2 = 1e30;
const Star *minstar = NULL;
double dist2;
for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Star *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;}
   }
 return (Star *) minstar; // violates constness
}
template<class Star> Star* StarList<Star>::FindAndRemoveClosest(double X, double Y) 
{
double min_dist2 = 1e30;
const Star *minstar = NULL;
double dist2;
 StarIterator si_res;
for (StarIterator si = this->begin(); si!= this->end(); ++si) 
   { 
   const Star *s = *si;
   dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
   if (dist2 < min_dist2) { min_dist2 = dist2; minstar = s;
   si_res = si ;}
   }
 if ( minstar == NULL)
 return (Star *) minstar; 
 else
   {
     Star * star_res = new Star(*minstar);
     this->erase(si_res);
     return (Star *) star_res; // violates constness
   }
}

template<class Star>bool StarList<Star>::HasCloseNeighbor(double X, double Y, double maxdist, double mindist, double minflux) const
{
  double dist2;
  double mindist2 = mindist*mindist;
  double maxdist2 = maxdist*maxdist;
  for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Star *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      bool okflux = true ;
      if (minflux > 0 )
	okflux = (s->flux > minflux);
      if (okflux && (dist2 > mindist2) && (dist2 < maxdist2)) {
	//cerr << "HasCloseNeighbor at dist : " << sqrt(dist2) << " Neighbor: " ;
	//s->dump();
	return true;}
    }
  return false;
}  

template<class Star> Star* StarList<Star>::ClosestNeighbor(double X, double Y, double mindist) const
{
  const Star *minstar = NULL;
  double dist2;
  double mindist2 = mindist*mindist;
  double min_dist2 = 1e30;
  for (StarCIterator si = this->begin(); si!= this->end(); ++si) 
    { 
      const Star *s = (*si); 
      dist2 = (X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y);
      if ((dist2 > mindist2) && (dist2 < min_dist2)) { min_dist2 = dist2; minstar = s;}
    }
  return (Star *) minstar;  // violates constness
}


template<class Star> int StarList<Star>::NumberOfNeighbors(const double &X, const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;

  for (StarCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Star *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  ++nstars;
	}
    }
  return nstars;
}

template<class Star> int StarList<Star>::AllNeighbors(StarList &NeighborList, const double &X, 
						      const double &Y, const double &distmax) const
{
  int nstars = 0;
  double dist2 = distmax*distmax;
  NeighborList.ClearList();
  for (StarCIterator it = this->begin(); it!= this->end(); ++it) 
    { 
      const Star *s  = *it; 
      if ((X - s->x)*(X - s->x) + (Y - s->y)*(Y - s->y) < dist2) 
	{
	  NeighborList.push_back(new Star(*s));
	  ++nstars;
	}
    }
  return nstars;
}

template<class Star>
bool DecreasingFlux(const Star *S1, const Star *S2)
{
  return (S1->flux > S2->flux);
}

template<class Star>void StarList<Star>::FluxSort()
{
  this->sort(&DecreasingFlux<Star>);
}

template<class Star>void StarList<Star>::CutTail(const int NKeep)
{
int count = 0;
StarIterator si;
for (si = this->begin(); si != this->end() && count < NKeep; ++count, ++si);
while ( si != this->end() ) {si = this->erase(si);}
}


template<class Star>void StarList<Star>::ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const
{
  for (StarCIterator s= this->begin(); s!= this->end(); ++s)
    {
      const Star *st  = *s;
      if (aFrame.InFrame(*st))
	{
	  Star *copy = new Star(*st);
	  Out.push_back(copy);
	}
    }
}

template<class Star>void StarList<Star>::CutEdges(const Frame &aFrame, float mindist) 
{
  for (StarIterator si= this->begin(); si!= this->end();)
    {
      if (aFrame.MinDistToEdges(**si) < mindist)
	{
	  si = this->erase(si); 
	}
      else
	si++;
    }
}

template<class Star>void StarList<Star>::CopyTo(StarList<Star> &Copy) const
{
  Copy.ClearList();
  Copy.GlobVal() = this->GlobVal();
  StarCIterator si;
  for (si = this->begin(); si != this->end(); ++si) Copy.push_back(new Star(*(*si)));
}


}}} // end of namespaces

#endif /* STARLIST__CC */
