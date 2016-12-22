#ifndef CCDIMAGE__H
#define CCDIMAGE__H

#include <string>
#include <list>
#include <vector>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
//#include "lsst/jointcal/Jointcal.h"


namespace lsst {
namespace jointcal {

typedef int VisitIdType;
typedef int CcdIdType;

/**
 * Handler of an actual image from a single CCD.
 * NOTE: could possibly be replaced with a subclass of afw.image.Exposure?
 */
class CcdImage : public RefCount
{
private:

    Frame imageFrame; // in pixels

    MeasuredStarList wholeCatalog; // the catalog of measured objets
    MeasuredStarList catalogForFit;

    CountedRef<BaseTanWcs> readWcs; // i.e. from pix to sky
    CountedRef<Gtransfo> inverseReadWcs; // i.e. from sky to pix

    // The following ones should probably be mostly removed.
    CountedRef<Gtransfo> CTP2TP; // go from CommonTangentPlane to this tangent plane.
    CountedRef<Gtransfo> TP2CTP; // reverse one
    CountedRef<Gtransfo> pix2CommonTangentPlane;// pixels -> CTP
    CountedRef<Gtransfo> pix2TP;

    CountedRef<Gtransfo> sky2TP;

    std::string riName;
    std::string riDir;
    std::string instrument;
    CcdIdType _ccdId;
    VisitIdType _visit;

    double airMass; // airmass value.
    double fluxCoeff; // coefficient to convert ADUs to ADUs/sec at airmass 1
    double mjd; // modified julian date
    double toadsZeroPoint;
    double elixirZP;
    double photk;
    double photc;
    double zp;
    double psfzp;
    std::string dateObs;
    // refraction
    double sineta, coseta, tgz, hourAngle; // eta : parallactic angle, z: zenithal angle (X = 1/cos(z))

    std::string _filter;
    std::string snlsgrid; // our grid corrections
    std::string flatcvmap; // a multiplicative map to apply to the fluxes
    int    index;
    int    expindex;

    Point  commonTangentPoint;

    void LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat, const std::string &fluxField);

public:

    CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
             const Point &CommonTangentPoint,
             const PTR(lsst::afw::image::TanWcs) wcs,
             const PTR(lsst::afw::image::VisitInfo) visitInfo,
             const lsst::afw::geom::Box2I &bbox,
             const std::string &filter,
             const PTR(lsst::afw::image::Calib) calib,
             const int &visit,
             const int &ccd,
             const std::string &fluxField );
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
    { return pix2CommonTangentPlane.get();}

    //!
    const Gtransfo* CommonTangentPlane2TP() const
    { return CTP2TP.get();}

    //!
    const Gtransfo* TP2CommonTangentPlane() const
    { return TP2CTP.get();}

    //!
    const Gtransfo* Pix2TangentPlane() const
    { return pix2TP.get();}

    //!
    const Gtransfo* Sky2TP() const
    { return sky2TP.get();}

    //! returns ccd ID
    int getCcdId() const { return _ccdId;}

    //! instrument (TOADINST fits pseudo-key)
    std::string Instrument() const {return instrument;}

    //! returns visit ID
    VisitIdType getVisit() const { return _visit;}

    //!  Airmass
    double getAirMass() const {return airMass;}

    //! Date Obs
    std::string getDateObs() const { return dateObs; }

    //! Julian Date
    double getMjd() const { return mjd; }

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
    Point RefractionVector() const {return Point(tgz*coseta, tgz*sineta);}

    //!conversion from ADU to ADU/sec at airmass=1
    double FluxCoeff() const { return fluxCoeff;}

    //! return the CcdImage filter name
    std::string getFilter() const { return _filter;}

    //! SNLS grid
    std::string SNLSGrid() const { return snlsgrid; }

    //! correction map to convert from one set of fluxes to another
    std::string FlatCVMap() const { return flatcvmap; }

    //void SetPix2TangentPlane(const Gtransfo *);

    //! the wcs read in the header. NOT updated when fitting.
    const Gtransfo *ReadWCS() const {return readWcs.get();}

    //! the inverse of the one above.
    const Gtransfo *InverseReadWCS() const {return inverseReadWcs.get();}

    //! Frame in pixels
    const Frame& ImageFrame() const { return imageFrame;}

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


/**
 * A list of CcdImage. Usually produced by Associations.
 */
//class CcdImageList : public std::list<CountedRef<CcdImage> >
class CcdImageList : public std::list<std::shared_ptr<CcdImage> >
{
public:

    template<class Accept> CcdImageList SubList(const Accept &OP) const
    {
        CcdImageList out;
        for (const_iterator i = begin(); i != end() ; ++i)
            if (OP(**i)) out.push_back(*i);
        return out;
    }
};


typedef CcdImageList::iterator CcdImageIterator;
typedef CcdImageList::const_iterator CcdImageCIterator;

}
} // end of namespaces

#endif /* CCDIMAGE__H */
