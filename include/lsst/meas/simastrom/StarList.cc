#ifndef STARLIST__CC
#define STARLIST__CC
#include <iostream>
#include <iomanip>
#include <istream>
#include <string>
#include <cstdio>
#include <string.h> // for strstr


#include "lsst/meas/simastrom/StarList.h"
#include "lsst/meas/simastrom/Frame.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace meas {
namespace simastrom {



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
	  std::cout << "WARNING: ignoring @lines " << std::endl;
	  //	  glob.ProcessLine(buff);
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
	  this->push_back(s); 
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
  std::ifstream rd(FileName.c_str());
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
  pr  << std::resetiosflags(std::ios::scientific) ;
  // pr  << setiosflags(0) ;
  pr  << std::setiosflags(std::ios::fixed) ;
  int oldprec = pr.precision();
  pr<< std::setprecision(8);

#if (0)
  // write GlobalValues if any
  std::vector<std::string> globs = glob.OutputLines();
  for (unsigned k = 0; k < globs.size(); ++k)
    pr << '@' << globs[k] << std::endl;
#endif

  // cannot use front() to detect emptyness
  if (this->begin() == this->end()) // empty std::list, and faster than (size() == 0)
    {
      Star dummy;
      dummy.WriteHeader(pr);
    }
  else this->front()->WriteHeader(pr);
  for (auto it= this->begin(); it!=this->end() ; it++ )
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
for (auto s= this->begin(); s!= this->end() && (count < NHead); ++s, count++)
  {
  Star *copy = new Star(*(*s));  /* to be changed if switch to Star rather than pointers to Stars */
  Out.push_back(copy);
  }
}

#if 0
template<class Star>
bool DecreasingFlux(const Star *S1, const Star *S2)
{
  return (S1->flux > S2->flux);
}
#endif

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
  for (auto s= this->begin(); s!= this->end(); ++s)
    {
      auto &st  = *s;
      if (aFrame.InFrame(*st))
	{
	  Star *copy = new Star(*st);
	  Out.push_back(copy);
	}
    }
}

template<class Star>void StarList<Star>::CutEdges(const Frame &aFrame, float mindist) 
{
  for (auto si= this->begin(); si!= this->end();)
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
  //  Copy.GlobVal() = this->GlobVal();
  for (auto si = this->begin(); si != this->end(); ++si) 
    Copy.push_back(new Star(*(*si)));
}


}}} // end of namespaces

#endif /* STARLIST__CC */
