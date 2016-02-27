// -*- C++ -*-
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "lsst/jointcal/RefStar.h"

namespace lsst {
namespace jointcal {


RefStar::RefStar(const BaseStar &B, const Point &RaDec) 
  : BaseStar(B), index(0)
{
    
  fittedStar = (FittedStar*) NULL;
  raDec = RaDec;
}

void RefStar::SetFittedStar(FittedStar &F)
{ 
  fittedStar = &F;
  if (fittedStar) fittedStar->SetRefStar(this);
}


double RefStar::Flux(int band) const
{
  if(band<0 && band>=10) return FP_INFINITE;
  return refFlux[band];
}


void RefStar::AssignRefFluxes(std::vector<double> const& reffluxes)
{
  refFlux.clear();
  copy(reffluxes.begin(), reffluxes.end(), back_inserter(refFlux));
}


//#include <starlist.cc>
/** RefStarList ***/
//template class StarList<RefStar>; /* to force instanciation */



BaseStarList& Ref2Base(RefStarList &This)
{
  return (BaseStarList&) This;
}
  
BaseStarList* Ref2Base(RefStarList *This)
{
  return (BaseStarList*) This;
}
  
const BaseStarList& Ref2Base(const RefStarList &This)
{
  return (const BaseStarList &) This;
}
  
const BaseStarList* Ref2Base(const RefStarList *This)
{
  return (BaseStarList*) This;
}
  


/********* RefStarStuple *************/
RefStarTuple::RefStarTuple(const std::string &FileName) 
  : stream(FileName.c_str())
{
  stream << "# ra : ref RA" << std::endl
	 << "# dec : ref DEC" << std::endl
	 << "# xref: deg in TP (from catalog)" << std::endl
	 << "# yref : " << std::endl
	 << "# xfit: deg in TP (from fit)" << std::endl
	 << "# yfit : " << std::endl
	 << "# dx : xref-xfit" << std::endl
	 << "# dy : yref-yfit" << std::endl
	 << "# refmag : ref mag" << std::endl
	 << "# fmag : f mag" << std::endl
	 << "# nm : number of measurements" << std::endl
	 << "# end " << std::endl;
    ;
  
}

void RefStarTuple::AddEntry(const RefStar &R, const FittedStar &F)
{
  std::ios::fmtflags  old_flags =  stream.flags();
  stream << std::setprecision(10);
  stream << R.Ra() << ' '
	 << R.Dec() << ' '
	 << R.x << ' ' 
	 << R.y << ' ' 
	 << F.x << ' '
	 << F.y << ' '
	 << R.x - F.x << ' ' 
	 << R.y - F.y << ' ' 
    ;
  stream << std::setprecision(5);
  stream << R.flux << ' '
	 << F.Mag() << ' '
	 << F.MeasurementCount() << ' ' 
	 << std::endl;

  stream.flags(old_flags);
}  

}} // end of namespaces

