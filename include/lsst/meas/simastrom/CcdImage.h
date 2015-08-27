#ifndef CCDIMAGE__H
#define CCDIMAGE__H

#include <string>
#include <list>
#include <vector>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/meas/simastrom/MeasuredStar.h"
#include "lsst/meas/simastrom/Gtransfo.h"
#include "lsst/meas/simastrom/Frame.h"


namespace lsst {
namespace meas {
namespace simastrom {



typedef int ShootIdType;


void SetZpKey(const std::string &AKey);

//! handler of an actual image from a single CCD
/*! requires an in-depth cleanup */
class CcdImage : public RefCount
{
 private:

  Frame imageFrame; // in pixels
  // wholeCatalog is just store the catalog of selected sources 
  //  lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> wholeCatalog;
  
  MeasuredStarList wholeCatalog; // the catalog of measured objets
  MeasuredStarList catalogForFit;

  // these 2 transfos are NOT updated when fitting
//  Gtransfo *readWcs; // i.e. from pix to sky
//  Gtransfo *inverseReadWcs; // i.e. from sky to pix  
  BaseTanWcs *readWcs; // i.e. from pix to sky
  Gtransfo *inverseReadWcs; // i.e. from sky to pix

  // The following ones should probably be mostly removed.
  Gtransfo *CTP2TP; // go from CommonTangentPlane to this tangent plane.
  Gtransfo *TP2CTP; // reverse one
  Gtransfo *pix2CommonTangentPlane;// pixels -> CTP
  Gtransfo *pix2TP;
  
  Gtransfo* sky2TP;
  
  std::string riName;
  std::string riDir;
  std::string instrument;
  int chip; // CCD number
  ShootIdType shoot; // Same value for all CcdImages from the same exposure
  unsigned bandRank; // some incremental band indicator.
  

  double expTime; // seconds
  double airMass; // airmass value.
  double fluxCoeff; // coefficient to convert ADUs to ADUs/sec at airmass 1
  double jd; // julian date
  double toadsZeroPoint;
  double elixirZP;
  double photk;
  double photc;
  double zp;
  double psfzp;
  //  double seeing;
  //  double gfseeing;
  //  double sigmaback;
  std::string dateObs;
  // refraction
  double sineta,coseta,tgz,hourAngle;  // eta : parallactic angle, z: zenithal angle (X = 1/cos(z))
  
  std::string band;
  std::string flatName;
  //  std::string flat; // full flat name
  //  std::string baseflat; // full baseflat name
  std::string cfhtscatter; // full scatter name
  std::string snlsgrid; // our grid corrections
  std::string flatcvmap; // a multiplicative map to apply to the fluxes
  int    bandIndex;
  int    index;
  int    expindex;
  
  Point  commonTangentPoint;

  void LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat);

 public:

  CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri, 
    const Point &CommonTangentPoint, const PTR(lsst::afw::image::TanWcs) wcs, const PTR(lsst::daf::base::PropertySet) meta,
    const lsst::afw::geom::Box2I &bbox, const std::string &filter, const PTR(lsst::afw::image::Calib) calib,
    const int &visit, const int &ccd, const std::string &ccdImage );


    
#ifdef TO_BE_FIXED 
  //!
  CcdImage(const ReducedImage &Ri, const Point &CommonTangentPoint, const CatalogLoader * LoadIt);
#endif

  //!
  std::string Name() const { return riName;}
  
  //!
  std::string Dir() const { return riDir; }

  //!
  const MeasuredStarList &WholeCatalog() const { return wholeCatalog;}

  //!
  const MeasuredStarList & CatalogForFit() const { return catalogForFit;}

  //!
  MeasuredStarList & CatalogForFit()  { return catalogForFit;}

  //! 
  const Gtransfo* Pix2CommonTangentPlane() const 
    { return pix2CommonTangentPlane;}

  //! 
  const Gtransfo* CommonTangentPlane2TP() const 
    { return CTP2TP;}
  
  //! 
  const Gtransfo* TP2CommonTangentPlane() const 
    { return TP2CTP;}
  
  //! 
  const Gtransfo* Pix2TangentPlane() const 
    { return pix2TP;}

  //! 
  const Gtransfo* Sky2TP() const 
    { return sky2TP;}
  
  //! returns chip ID
  int Chip() const { return chip;}

  //! instrument (TOADINST fits pseudo-key)
  std::string Instrument() const {return instrument;}

  //! some incremental band rank. Is used to incrementally index bands in a sample of input images. Different from BandIndex() 
  unsigned BandRank() const {return bandRank;}

  //! returns seeing
  //  double Seeing() const { return seeing;}
  
  //! returns gfseeing
  //  double GFSeeing() const { return gfseeing;}
  
  //! returns sigma back
  //  double SigmaBack() const { return sigmaback;}
  
  //! returns shoot ID
  ShootIdType Shoot() const { return shoot;}
  
  //! Exposure time (s)
  double ExpTime() const { return expTime;}
  
  //!  Airmass
  double AirMass() const {return airMass;}
  
  //! Date Obs
  std::string DateObs() const { return dateObs; }

  //! Julian Date
  double JD() const { return jd; }

  
  //!Elixir ZP (applies to fluxes in ADU/sec at airmass 1).
  double ElixirZP() const { return elixirZP;}

  //!zp from the Fits key set by SetZpKey(std::string)
  double ZP() const { return zp;}

  //!zp from the psf zp file, returns 0 if not present
  double PSFZP() const { return psfzp;}
  
  //! absorption term
  double PhotK() const { return photk; }
  
  //! original ZP 
  double PhotC() const { return photc; }

  //! 
  double HourAngle() const { return hourAngle; }

  //! Parallactic angle
  double SinEta() const { return sineta; }

  //! Parallactic angle
  double CosEta() const { return coseta; }

  //! Parallactic angle
  double TanZ() const { return tgz; }

  //!
  Point ParallacticVector() const {return Point(tgz*coseta, tgz*sineta);}

  //!conversion from ADU to ADU/sec at airmass=1
  double FluxCoeff() const { return fluxCoeff;}
  
  //! return the CcdImage band name
  std::string Band() const { return band;}
  
  //! return the CcdImage band index. This is a static index that mostly turns a letter (e.g. 'g') into a number (e.g. 2). Different from BandRank() 
  int BandIndex() const { return bandIndex; }
  
  //! Flat used to flatfield
  std::string FlatName() const { return flatName;}

  //! Full path of the scatter corrections
  std::string CFHTScatter() const { return cfhtscatter; }
  
  //! SNLS grid
  std::string SNLSGrid() const { return snlsgrid; }
  
  //! correction map to convert from one set of fluxes to another
  std::string FlatCVMap() const { return flatcvmap; }
  
  void SetPix2TangentPlane(const Gtransfo *);
  
  //! the wcs read in the header. NOT updated when fitting.
  const Gtransfo *ReadWCS() const {return readWcs;}
  
  //! the inverse of the one above.
  const Gtransfo *InverseReadWCS() const {return inverseReadWcs;}
  
  //! Frame in pixels
  const Frame& ImageFrame() const { return imageFrame;}
  
  //! Frame on sky
  Frame RaDecFrame() const;
  
  //! Fitted Ccd object (contain the refscale parameters)
  //  void             SetFittedCcd(FittedCcd* ccd) { if(ccd) fittedccd=ccd; }
  //  FittedCcd*       GetFittedCcd() { return fittedccd; }
  //  FittedCcd const* GetFittedCcd() const { return fittedccd; }
  
  //! returns wether "this" overlaps with Other.
  //  bool Overlaps(const CcdImage &Other) const;
  
  //! CcdImage index
  int     Index() const { return index; }
  void    SetIndex(int idx) { index = idx; }
  
  //! Exposure Index
  int     ExpIndex() const { return expindex; }
  void    SetExpIndex(int idx) { expindex = idx; }
  
  //! Common Tangent Point
  Point const&       CommonTangentPoint() const { return commonTangentPoint; }
  
 private:
  CcdImage(const CcdImage &); // forbid copies


};


/********* CcdImageList *************/


class CcdImageList : public std::list<CountedRef<CcdImage> >
{
  public:
  
  //! 
  std::list<std::string> DateObs() const;
  
  //! 
  std::list<std::string> Bands() const;
  
  //! 
  double       MeanAirmass() const;
  
  //! 
  template<class Accept> CcdImageList SubList(const Accept &OP) const
    {
      CcdImageList out;
      for (const_iterator i = begin(); i != end() ; ++i)
	if (OP(**i)) out.push_back(*i);
      return out;
    }
  
  // find the matching image. Chip==-1 means any chip
  double AirMass(const int Shoot, const int Chip = -1) const;
};


typedef CcdImageList::iterator CcdImageIterator;
typedef CcdImageList::const_iterator CcdImageCIterator;

}}} // end of namespaces

#endif /* CCDIMAGE__H */
