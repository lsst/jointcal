#ifndef CCDIMAGE__H
#define CCDIMAGE__H

#include <string>
#include <list>
#include <vector>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"

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

    MeasuredStarList _wholeCatalog; // the catalog of measured objets
    MeasuredStarList _catalogForFit;

    CountedRef<BaseTanWcs> readWcs; // i.e. from pix to sky
    CountedRef<Gtransfo> inverseReadWcs; // i.e. from sky to pix

    // The following ones should probably be mostly removed.
    CountedRef<Gtransfo> CTP2TP; // go from CommonTangentPlane to this tangent plane.
    CountedRef<Gtransfo> TP2CTP; // reverse one
    CountedRef<Gtransfo> pix2CommonTangentPlane;// pixels -> CTP
    CountedRef<Gtransfo> pix2TP;

    CountedRef<Gtransfo> sky2TP;

    std::string name;
    CcdIdType _ccdId;
    VisitIdType _visit;

    lsst::afw::coord::IcrsCoord boresightRaDec;
    double airMass; // airmass value.
    double fluxCoeff; // coefficient to convert ADUs to ADUs/sec at airmass 1
    double mjd; // modified julian date
    PTR(lsst::afw::image::Calib) _calib;
    std::string dateObs;
    // refraction
    double sineta, coseta, tgz, hourAngle; // eta : parallactic angle, z: zenithal angle (X = 1/cos(z))

    std::string _filter;

    Point  _commonTangentPoint;

    void LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat, const std::string &fluxField);

public:

    CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Ri,
             const PTR(lsst::afw::image::TanWcs) wcs,
             const PTR(lsst::afw::image::VisitInfo) visitInfo,
             const lsst::afw::geom::Box2I &bbox,
             const std::string &filter,
             const PTR(lsst::afw::image::Calib) calib,
             const int &visit,
             const int &ccd,
             const std::string &fluxField );

    //! Return the name that identifies this ccdImage.
    std::string getName() const { return name;}

    /**
     * @brief      Gets the as-read catalog.
     *
     * @return     The whole catalog.
     */
    const MeasuredStarList &getWholeCatalog() const { return _wholeCatalog;}

    //@{
    /**
     * @brief      Gets the catalog to be used for fitting, which may have been cleaned-up.
     *
     * @return     The catalog for fitting.
     */
    const MeasuredStarList & getCatalogForFit() const { return _catalogForFit;}
    MeasuredStarList & getCatalogForFit()  { return _catalogForFit;}
    //@}

    /**
     * @brief      Sets the common tangent point and computes necessary transforms.
     *
     * @param[in]  commonTangentPoint  The common tangent point of all ccdImages (decimal degrees).
     */
    void setCommonTangentPoint(const Point &commonTangentPoint);

    /**
     * @brief      Gets the common tangent point, shared between all ccdImages.
     *
     * @return     The common tangent point of all ccdImages (decimal degrees).
     */
    Point const& getCommonTangentPoint() const { return _commonTangentPoint; }

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

    //! returns visit ID
    VisitIdType getVisit() const { return _visit;}

    //!  Airmass
    double getAirMass() const {return airMass;}

    //! Date Obs
    std::string getDateObs() const { return dateObs; }

    //! Julian Date
    double getMjd() const { return mjd; }

    //!Return the exposure's photometric calibration
    PTR(lsst::afw::image::Calib) getCalib() { return _calib; }

    /**
     * @brief      Gets the boresight RA/Dec.
     */
    lsst::afw::coord::IcrsCoord getBoresightRaDec() { return boresightRaDec; }

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

    //! the wcs read in the header. NOT updated when fitting.
    const Gtransfo *ReadWCS() const {return readWcs.get();}

    //! the inverse of the one above.
    const Gtransfo *InverseReadWCS() const {return inverseReadWcs.get();}

    //! Frame in pixels
    const Frame& ImageFrame() const { return imageFrame;}

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
        for (auto const &i: *this)
            if (OP(*i)) out.push_back(i);
        return out;
    }

    /**
     * @brief      Return the number of ccdImages with a populated catalogForFit list.
     */
    int sizeValidForFit() const
    {
        return std::count_if(this->begin(), this->end(),
                             [](std::shared_ptr<CcdImage> const &item)
                             {return item->getCatalogForFit().size() > 0;});
    }
};


typedef CcdImageList::iterator CcdImageIterator;
typedef CcdImageList::const_iterator CcdImageCIterator;

}
} // end of namespaces

#endif /* CCDIMAGE__H */
