// -*- C++ -*-
#include <cmath>
#include <vector>

#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/StarList.cc"
#include "lsst/jointcal/CcdImage.h"

//#include "preferences.h"
//#include "ccdimage.h"
#include "assert.h" // for assert


namespace lsst {
namespace jointcal {

  /* Interesting fields of the stack catalogs :
 'base_SdssCentroid_x'
 'base_SdssCentroid_y'
 'base_SdssCentroid_xSigma'
 'base_SdssCentroid_ySigma'

  We miss the xy uncertainty term.
  We can cook it up from the sdss shape:
  'base_SdssShape_xx'
  'base_SdssShape_yy'
  'base_SdssShape_xy'

  for fluxes, we might use :
  'base_CircularApertureFlux_2_flux'
  'base_CircularApertureFlux_2_fluxSigma'

   where the '2' should be read from the environment.
  */






//#define WRITE_STARS

#ifdef TO_BE_FIXED
MeasuredStar::MeasuredStar( const SEStar &S) :
  BaseStar(S),
  mag(0.), wmag(0.), eflux(0.), aperrad(0.), ccdImage(0),
  fittedStar((const FittedStar *)NULL),
  valid(true)
{
  eflux = S.EFlux();
  if (vx == 1 ) // means no position errors in input catalogs
    {
      cout << "[ INFO : missing position uncertainties in input data] " << endl;
      double sig = S.Fwhm()/2.36*sqrt(eflux/flux)/5.; // totally ad'hoc expression
      vx = vy = sig*sig;
      vxy = 0;
    }
  mag = -1;
  //  InitialStarRef = &S;
  
  // Some (temporary) User values
  usrVals.push_back(S.Mxx());
  usrVals.push_back(S.Myy());
  usrVals.push_back(S.Mxy());
  usrVals.push_back(S.A());
  usrVals.push_back(S.B());
  usrVals.push_back(S.Gyr_Angle());
}
#endif


BaseStarList& Measured2Base(MeasuredStarList &This)
{
  return (BaseStarList&) This;
}

BaseStarList* Measured2Base(MeasuredStarList *This)
{
  return (BaseStarList*) This;
}

const BaseStarList& Measured2Base(const MeasuredStarList &This)
{
  return (const BaseStarList &) This;
}

const BaseStarList* Measured2Base(const MeasuredStarList *This)
{
  return (BaseStarList*) This;
}


//! StarList ascii IO's

  std::string MeasuredStar::WriteHeader_(std::ostream & pr ,
					 const char* i) const
{
  //  string format = BaseStar::WriteHeader_(pr,i);
  std::string format = BaseStar::WriteHeader_(pr,i);
  if(i==NULL) i = "";
  
  pr << "# eflux" << i << " : " << std::endl
     << "# mag" << i << " : " << std::endl
     << "# wmag" << i << " : " << std::endl
     << "# ccd"<<i<<" : chip number" << std::endl
     << "# expo"<<i<<" : exposure number" << std::endl
     << "# airmass"<<i<<" :  "<<std::endl
     << "# band"<<i<< " : ugriz=01234 " << std::endl
     << "# valid" << i << " : 0=outlier 1=ok ?" << std::endl;
  
  format += " MeasuredStar 1 ";
  //  if(i && InitialStarRef)
  //    {
  //      string i2 = "1"+string(i);
  //      format += InitialStarRef->WriteHeader_(pr,i2.c_str());
  //    }
  return format;
}


BaseStar* MeasuredStar::read(std::istream& s, const char* format)
{
  MeasuredStar* ret = new MeasuredStar();
  //  ret->BaseStar::read_it(s, format);
  ret->BaseStar::read_it(s, format);
  double tmp;
  s >> ret->eflux
    >> ret->mag
    >> ret->wmag
    >> tmp // chip
    >> tmp // shoot
    >> tmp // airmass
    >> tmp // band
    >> ret->valid;
  return ret;
}


void MeasuredStar::writen(std::ostream& s) const
{
  //  static map<string,int> bandNumber;
  //  static bool called = false;
  //  if (!called)
  //    {
  //      called = true;
  //      bandNumber["u"]=0;
  //      bandNumber["g"]=1;
  //      bandNumber["r"]=2;
  //      bandNumber["i"]=3;
  //      bandNumber["z"]=4;
  //    }
  
  BaseStar::writen(s);

  s << eflux << ' '
    << (std::isnan(mag)? -1 : mag) << ' '
    << wmag << ' ';
  
  if(ccdImage)
    std::cout << "ccdImage ! " << ccdImage << std::endl;
  
  if(ccdImage)
    s << ccdImage->Chip() << ' '
      << ccdImage->Shoot() << ' '
      << ccdImage->AirMass() << ' '
      << ccdImage->BandIndex() << ' ';
  else
    s << 0 << ' '
      << 0 << ' '
      << 0 << ' '
      << -1 << ' ';
  
  s << valid << ' ';

}



/******* MeasuredStarList *********/





void MeasuredStarList::SetZeroPoint(const double &ZP)
{
  zeroPoint = ZP;
  for (MeasuredStarIterator i= begin(); i != end(); ++i)
    {
      MeasuredStar &m = **i;
      m.mag = -2.5*log10(m.flux)+zeroPoint;
    }
}


void  MeasuredStarList::SetCcdImage(const CcdImage *C)
{
  for (MeasuredStarIterator i= begin(); i != end(); ++i)
      (*i)->SetCcdImage(C);
}
  




template class StarList<MeasuredStar>; // to force instanciation


#ifdef WE_WILL_USE_SOME_OTHER_MECHANISM

/************* routines related to loading input catalogs *****/


/*
   did not focus on flexibility... The trivial way to change
   cuts on loading scheme is to provide a home made catalog loader
   to MeasuredStarList constructor.  Associations' constructor
   propagates it all the way through.
*/


/* selection routine for SExtractor input catalogs */
bool KeepSEStar(const SEStar &S)
{
  return (S.Flag() ==0
	  && S.FlagBad() == 0
	  &&  S.EFlux() > 0
          && S.Fwhm() > 0
	  && S.flux > Preferences().minSigToNoise*S.EFlux());
}


bool KeepAperSEStar(const AperSEStar &S)
{
 if (S.Flag() == 0
     && S.FlagBad() == 0
     && S.neighborContamination == false)
   {
     //     const AperSEStar::Aperture &aper = S.apers[Preferences().aperRank];
     const Aperture &aper = S.apers[Preferences().aperRank];
     return (aper.eflux > 0 &&
	     aper.flux > aper.eflux * Preferences().minSigToNoise);
   }
 return false;
}




bool SECatalogLoader_(const ReducedImage &Ri, const CcdImage &ccdImage,
			    MeasuredStarList &List)
{
  if (!&ccdImage) {} // warning killer
  List.clear();
  SEStarList seList(Ri.CatalogName());
  for (SEStarCIterator i= seList.begin(); i != seList.end(); ++i)
    {
      const SEStar &se = **i;
      if (KeepSEStar(se))
	List.push_back(new MeasuredStar(se));
    }
  List.SetZeroPoint(Ri.ZeroPoint());
  return true;
}



bool AperSECatalogLoader_(const ReducedImage &Ri, const CcdImage &ccdImage,
				MeasuredStarList &List)
{
  if (!&ccdImage) {} // warning killer
  List.clear();
  string fileName(Ri.AperCatalogName());
  AperSEStarList al(fileName);
  unsigned rank = Preferences().aperRank;
  //  double seeing = al.GlobVal().Value("SEEING");
  double seeing = al.GlobVal().getDoubleValue("SEEING");
  //  vector<double> rads = al.GlobVal().Values("RADIUS");
  vector<double> rads = al.GlobVal().getDoubleValues("RADIUS");
  if (rads.size() < rank + 1)
    {
      cerr << " not enough apertures in " << fileName << std::endl;
      return false;
    }
  cout << " reading aper fluxes in " << fileName << " for aperture " <<
    rads[rank] << " pixels = " << rads[rank]/seeing << " seeing(r.m.s)" << std::endl;
  for (AperSEStarCIterator i = al.begin(); i != al.end(); ++i)
    {
      const AperSEStar &a = **i;
      if (a.apers[rank].nbad != 0) continue;
      if (KeepAperSEStar(a))
	{
	  MeasuredStar *m = new MeasuredStar(a);
	  m->flux = a.apers[rank].flux;
	  m->eflux = a.apers[rank].eflux;
	  List.push_back(m);
	}
    }
  List.SetZeroPoint(Ri.ZeroPoint());
  return true;
}



CatalogLoader * CatalogLoader::getDefaultCatalogLoader(){
  switch (Preferences().catalogToRead)
    {
    case SE_CATALOG :
      return new CatalogLoader(&SECatalogLoader_);
      break;
    case APER_CATALOG :
      return new CatalogLoader(&AperSECatalogLoader_);
      break;
    default :
    std::cerr << " DONT KNOW THIS kind of catalog "
	   << Preferences().catalogToRead << std::endl
	   << " allowed values in Preferences().catalogToRead "
	   << " and CATALOG_TO_READ key in datacards " << std::endl
	   << APER_CATALOG << ',' << SE_CATALOG << std::endl;
      throw(PolokaException(" Unkown catalog type in CatalogLoader::getDefaultCatalogLoader "));
    }
  return NULL;
}

#endif /* WE_WILL_USE_SOME_OTHER_MECHANISM */

}} // end of namespaces






