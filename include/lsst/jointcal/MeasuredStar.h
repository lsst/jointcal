// -*- C++ -*-
#ifndef MEASUREDSTAR__H
#define MEASUREDSTAR__H

#include <iostream>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {



class CcdImage;

/*! \file */

//! objects measured on actual images. Coordinates and uncertainties are expressed in pixel image frame. Flux expressed in ADU/s.
class MeasuredStar : public BaseStar
{
  public :
  double mag;
  double wmag;
  double eflux;
  double aperrad;
  double chi2;


  const CcdImage *ccdImage;
  std::vector<double> usrVals;

  private :

    CountedRef<const FittedStar> fittedStar;
    bool   valid;



  public :

    //!
  MeasuredStar()
    : BaseStar(),
      mag(0.), wmag(0.), eflux(0.), aperrad(0.),
      ccdImage(0),
      valid(true) {}

  MeasuredStar(const BaseStar &B, const FittedStar *F = nullptr) :
    BaseStar(B),
    mag(0.), wmag(0.), eflux(0.), aperrad(0.),
    ccdImage(0),
    valid(true)

  {
    //    InitialStarRef = &B;
  }

    //  MeasuredStar( const SEStar &S);

  void SetFittedStar(FittedStar *F)
      { if (F)  F->MeasurementCount()++; fittedStar = F;
      }

  double FluxSig() const { return eflux;}

  double Mag() const { return mag;}
  double AperRad() const { return aperrad; }

  //! the inverse of the mag variance
  double MagWeight() const { return (flux*flux/(eflux*eflux));}

  const FittedStar* GetFittedStar() const { return fittedStar.get();};

  const CcdImage &GetCcdImage()  const { return *ccdImage;};

  void setCcdImage(const CcdImage *C) { ccdImage = C;};

  //! Fits may use that to discard outliers
  bool  IsValid() const { return valid; }
  //! Fits may use that to discard outliers
  void  SetValid(bool v) { valid=v; }

  // No longer decrement counter of associated fitted star in destructor (P. El-Hage le 10/04/2012)
  // ~MeasuredStar() { if (fittedStar) fittedStar->MeasurementCount()--;}

  std::string WriteHeader_(std::ostream & pr = std::cout, const char* i = nullptr) const;
  void writen(std::ostream& s) const;
static BaseStar*  read(std::istream &s, const char* format);
};


typedef CountedRef<MeasuredStar> MeasuredStarRef;




/****** MeasuredStarList */


//! A list of MeasuredStar. They are usually filled in Associations::AddImage
class MeasuredStarList : public StarList<MeasuredStar> {

  double zeroPoint;

  public :
    MeasuredStarList() {};


  void SetZeroPoint(double ZP);

  void setCcdImage(const CcdImage *C);
};



typedef MeasuredStarList::const_iterator MeasuredStarCIterator;
typedef MeasuredStarList::iterator MeasuredStarIterator;
typedef CountedRef<MeasuredStar> MeasuredStarRef;

BaseStarList& Measured2Base(MeasuredStarList &This);
BaseStarList* Measured2Base(MeasuredStarList *This);
const BaseStarList& Measured2Base(const MeasuredStarList &This);
const BaseStarList* Measured2Base(const MeasuredStarList *This);


#ifdef TO_BE_REPLACED

typedef bool (CatalogLoader_(const ReducedImage &R, const CcdImage &, MeasuredStarList &List));

class CatalogLoader{
  CatalogLoader_ * loader;
public:
  CatalogLoader(CatalogLoader_  *l){
    loader = l;
  }
  bool operator()(const ReducedImage & ri,
			  const CcdImage & ccdim,
			  MeasuredStarList & stl) const{
    return load(ri, ccdim, stl);
  }
  virtual bool load(const ReducedImage & ri,
		    const CcdImage & ccdim,
		    MeasuredStarList & stl) const{
    return (*loader)(ri, ccdim, stl);
  }

  virtual ~CatalogLoader(){};

  static CatalogLoader * DefaultCatalogLoader;
  static CatalogLoader * getDefaultCatalogLoader();

protected:
  CatalogLoader(){}
};

//! routine that checks the quality of SEStar measurements and applies the signal2noise cut from datacards.
bool KeepSEStar(const SEStar &S);

//! select stars and normalizes fluxes (ADU/sec @ X=1)
extern const CatalogLoader &AperSELoader_Minimal_NormalizeFluxes;
extern const CatalogLoader &AperSELoader_SelectStars_NormalizeFluxes;
extern const CatalogLoader &SELoader_SelectStars_NormalizeFluxes;
extern const CatalogLoader &AperSELoader_NormalizeFluxes;
extern const CatalogLoader &AperSELoader_SelectStars_NormalizeFluxes_nosncut;
extern const CatalogLoader &CachedMeasuredStarLoader;
extern const CatalogLoader &MCMeasuredStarLoader;


void SelectFluxForFitCatalog(std::string const& filename);
//void ConvertCatalogFluxesToMagnitudes(bool convert /*, double zp=0. */);

#endif /* TO_BE_REPLACED */

}}  // end of namespaces

#endif /* MEASUREDSTAR__H */

