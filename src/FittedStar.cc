#include <iostream>
#include <iomanip>

#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/StarList.cc"


namespace lsst {
namespace jointcal {


// cannot be in fittedstar.h, because of "crossed includes"
FittedStar::FittedStar(const MeasuredStar &M) :
  BaseStar(M), mag(M.Mag()), emag(-1), col(0.), gen(-1), wmag(M.MagWeight()),
  indexInMatrix(-1), measurementCount(0), refStar(NULL),
  flux2(-1), fluxErr2(-1)
{
  fluxErr = M.eflux;
}


void FittedStar::SetRefStar(const RefStar *R)
{
  if (refStar != NULL && (R)) // Exception ??
    std::cerr << " FittedStar : " << *this
	      << " is already matched to an other RefStar " << std::endl
	      << " Clean up your lists " << std::endl;
  else refStar = R;
}

static double sq(const double &x) {return x*x;}


void FittedStar::AddMagMeasurement(const double &MagValue,
				   const double &MagWeight)
{
  mag = (mag*wmag+MagValue*MagWeight)/(wmag+MagWeight);
  wmag += MagWeight;
}





/************* FittedStarList ************************/



#ifdef DO_WE_NEED_IT
/* I am not sure that reading using the DicStar mechanism is a good idea.
   If we need persistence of FittedStarLists, we'll devise the I/O's.
   Pierre Astier (July 15)
*/
#include "dicstar.h"
//! read a list from a previous run
FittedStarList::FittedStarList(const std::string &FileName)
{
  DicStarList dl(FileName);
  for (DicStarIterator i = dl.begin(); i != dl.end(); ++i)
    {
      DicStar &ds = **i;
      // The file contains in principle a mag, and BaseStar expects
      FittedStar *fs = new FittedStar(ds);
      fs->SetMag(ds.flux);
      push_back(fs);
    }
}

#endif /* DO_WE_NEED_IT */
  



BaseStarList& Fitted2Base(FittedStarList &This)
{
  return (BaseStarList&) This;
}
 
BaseStarList* Fitted2Base(FittedStarList *This)
{
  return (BaseStarList*) This;
}
 
const BaseStarList& Fitted2Base(const FittedStarList &This)
{
  return (const BaseStarList &) This;
}
 
const BaseStarList* Fitted2Base(const FittedStarList *This)
{
  return (BaseStarList*) This;
}



void FittedStarList::WriteTuple(const std::string &FileName,
			  const Gtransfo &T,
			  const bool OnlyGoodStars)
{
  FittedStarTuple tuple(FileName);
  for (FittedStarCIterator i = begin(); i != end(); ++i)
    {
      const FittedStar &f = **i;
      if (OnlyGoodStars && f.flux < 0) continue;
      Point raDec = T.apply(f);
      tuple.AddEntry(f, raDec);
    }
}


/****************/
 
FittedStarTuple::FittedStarTuple(const std::string &FileName)
  : stream(FileName.c_str())
{
  stream << "# ra : " << std::endl
	 << "# dec : " << std::endl
	 << "# mag : " << std::endl
	 << "# nm : number of measurements" << std::endl
	 << "# end " << std::endl;
    ;
  
}





void FittedStarTuple::AddEntry(const FittedStar &F, const Point &RaDec)
{
  std::ios::fmtflags  old_flags =  stream.flags();
  stream << std::setprecision(10);
  stream << RaDec.x << ' '
	 << RaDec.y << ' '
    ;
  stream << std::setprecision(5);
  stream << F.Mag() << ' '
	 << F.MeasurementCount() << ' '
	 << std::endl;

  stream.flags(old_flags);
}



template class StarList<FittedStar>;

std::string FittedStar::WriteHeader_(std::ostream& pr, const char* i) const
{
 std::string format = BaseStar::WriteHeader_(pr, i);
  if(i==NULL) i = "";
  
  pr << "# mag"   << i << " : fitted magnitude" << std::endl;
  pr << "# emag"  << i << " : error on the fitted magnitude" << std::endl;
  pr << "# col"   << i << " : star color" << std::endl;
  pr << "# gen"   << i << " : star generation (very grid specific)" << std::endl;
  pr << "# wmag"  << i << " : some other field" << std::endl;
  pr << "# index" << i << " : star index" << std::endl;
  pr << "# measc" << i << " : measurement count" << std::endl;
  pr << "# flux2" << i << " : flux in some other band ()" << std::endl;
  pr << "# fluxErr" << i << " : flux error" << std::endl;
  pr << "# fluxErr2" << i << " : flux2 error" << std::endl;
  
  format += "FittedStar 1";
  return format;
}



 void FittedStar::writen(std::ostream& s) const
{
  s.setf(std::ios::scientific);
  s << std::setprecision(12);
  
  BaseStar::writen(s);
  s << mag     << " "
    << emag    << " "
    << col     << " "
    << gen     << " "
    << wmag    << " "
    << indexInMatrix   << " "
    << measurementCount << " "
    << flux2   << " "
    << fluxErr << " "
    << fluxErr2;
}


void FittedStar::read_it(std::istream& s, const char* format)
{
  BaseStar::read_it(s, format);
  s >> mag
    >> emag
    >> col
    >> gen
    >> wmag
    >> indexInMatrix
    >> measurementCount
    >> flux2
    >> fluxErr
    >> fluxErr2;
}


BaseStar* FittedStar::read(std::istream& s, const char* format)
{
  FittedStar* ret = new FittedStar();
  ret->read_it(s, format);
  return ret;
}


}} // end of namespaces
