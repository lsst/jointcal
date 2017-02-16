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
  indexInMatrix(-1), measurementCount(0), _refStar(nullptr),
  flux2(-1), fluxErr2(-1)
{
  fluxErr = M.eflux;
}


void FittedStar::setRefStar(const RefStar *refStar)
{
  if ((_refStar != nullptr) && (refStar != nullptr)) // TODO: should we raise an Exception in this case?
    // TODO: This message should be log.warn()
    std::cerr << " FittedStar : " << *this
	      << " is already matched to another RefStar. Clean up your lists" << std::endl
          << "old: " << *_refStar << std::endl
          << "new: " << *refStar << std::endl;
  else _refStar = refStar;
}

void FittedStar::AddMagMeasurement(double MagValue,
				   double MagWeight)
{
  mag = (mag*wmag+MagValue*MagWeight)/(wmag+MagWeight);
  wmag += MagWeight;
}





/************* FittedStarList ************************/


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
  for (auto const &fittedStar: *this)
    {
      if (OnlyGoodStars && fittedStar->flux < 0) continue;
      Point raDec = T.apply(*fittedStar);
      tuple.AddEntry(*fittedStar, raDec);
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
  if(i==nullptr) i = "";

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
