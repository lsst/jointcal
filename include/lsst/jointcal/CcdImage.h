// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CCD_IMAGE_H
#define LSST_JOINTCAL_CCD_IMAGE_H

#include <list>
#include <string>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"

namespace lsst {
namespace jointcal {

typedef std::list<std::shared_ptr<CcdImage>> CcdImageList;

typedef int VisitIdType;
typedef int CcdIdType;

/**
 * Handler of an actual image from a single CCD.
 * NOTE: could possibly be replaced with a subclass of afw.image.Exposure?
 */
class CcdImage {
private:
    Frame _imageFrame;  // in pixels

    MeasuredStarList _wholeCatalog;  // the catalog of measured objets
    MeasuredStarList _catalogForFit;

    std::shared_ptr<BaseTanWcs> _readWcs;       // i.e. from pix to sky
    std::shared_ptr<Gtransfo> _inverseReadWcs;  // i.e. from sky to pix

    // The following ones should probably be mostly removed.
    std::shared_ptr<Gtransfo> _CTP2TP;                  // go from CommonTangentPlane to this tangent plane.
    std::shared_ptr<Gtransfo> _TP2CTP;                  // reverse one
    std::shared_ptr<Gtransfo> _pix2CommonTangentPlane;  // pixels -> CTP
    std::shared_ptr<Gtransfo> _pix2TP;

    std::shared_ptr<Gtransfo> _sky2TP;

    std::string _name;
    CcdIdType _ccdId;
    VisitIdType _visit;

    lsst::afw::coord::IcrsCoord _boresightRaDec;
    double _airMass;  // airmass value.
    double _mjd;      // modified julian date
    std::shared_ptr<afw::image::PhotoCalib> _photoCalib;
    // refraction
    // eta : parallactic angle, z: zenithal angle (X = 1/cos(z))
    double _sineta, _coseta, _tgz;
    // Local Sidereal Time and hour angle of observation
    double _lstObs, _hourAngle;

    std::string _filter;

    Point _commonTangentPoint;

    void LoadCatalog(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> const &Cat,
                     std::string const &fluxField);

public:
    CcdImage(afw::table::SourceCatalog &record, std::shared_ptr<lsst::afw::image::TanWcs> wcs,
             std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, afw::geom::Box2I const &bbox,
             std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib, int visit,
             int ccd, std::string const &fluxField);

    /// No move or copy: each CCD image is unique to that ccd+visit, and Associations holds all CcdImages.
    CcdImage(CcdImage const &) = delete;
    CcdImage(CcdImage &&) = delete;
    CcdImage &operator=(CcdImage const &) = delete;
    CcdImage &operator=(CcdImage &&) = delete;

    //! Return the _name that identifies this ccdImage.
    std::string getName() const { return _name; }

    /**
     * @brief      Gets the as-read catalog.
     *
     * @return     The whole catalog.
     */
    MeasuredStarList const &getWholeCatalog() const { return _wholeCatalog; }

    //@{
    /**
     * @brief      Gets the catalog to be used for fitting, which may have been cleaned-up.
     *
     * @return     The catalog for fitting.
     */
    MeasuredStarList const &getCatalogForFit() const { return _catalogForFit; }
    MeasuredStarList &getCatalogForFit() { return _catalogForFit; }
    //@}

    /**
     * @brief      Sets the common tangent point and computes necessary transforms.
     *
     * @param[in]  commonTangentPoint  The common tangent point of all ccdImages (decimal degrees).
     */
    void setCommonTangentPoint(Point const &commonTangentPoint);

    /**
     * @brief      Gets the common tangent point, shared between all ccdImages.
     *
     * @return     The common tangent point of all ccdImages (decimal degrees).
     */
    Point const &getCommonTangentPoint() const { return _commonTangentPoint; }

    //!
    Gtransfo const *getPix2CommonTangentPlane() const { return _pix2CommonTangentPlane.get(); }

    //!
    Gtransfo const *getCommonTangentPlane2TP() const { return _CTP2TP.get(); }

    //!
    Gtransfo const *getTP2CommonTangentPlane() const { return _TP2CTP.get(); }

    //!
    Gtransfo const *getPix2TangentPlane() const { return _pix2TP.get(); }

    //!
    Gtransfo const *getSky2TP() const { return _sky2TP.get(); }

    //! returns ccd ID
    int getCcdId() const { return _ccdId; }

    //! returns visit ID
    VisitIdType getVisit() const { return _visit; }

    //!  Airmass
    double getAirMass() const { return _airMass; }

    //! Julian Date
    double getMjd() const { return _mjd; }

    //! Return the exposure's photometric calibration
    std::shared_ptr<afw::image::PhotoCalib> getPhotoCalib() { return _photoCalib; }

    /**
     * @brief      Gets the boresight RA/Dec.
     */
    lsst::afw::coord::IcrsCoord getBoresightRaDec() { return _boresightRaDec; }

    //!
    double getHourAngle() const { return _hourAngle; }

    //! Parallactic angle
    double getSinEta() const { return _sineta; }

    //! Parallactic angle
    double getCosEta() const { return _coseta; }

    //! Parallactic angle
    double getTanZ() const { return _tgz; }

    //!
    Point getRefractionVector() const { return Point(_tgz * _coseta, _tgz * _sineta); }

    //! return the CcdImage filter name
    std::string getFilter() const { return _filter; }

    //! the wcs read in the header. NOT updated when fitting.
    Gtransfo const *readWCS() const { return _readWcs.get(); }

    //! the inverse of the one above.
    Gtransfo const *getInverseReadWCS() const { return _inverseReadWcs.get(); }

    //! Frame in pixels
    Frame const &getImageFrame() const { return _imageFrame; }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CCD_IMAGE_H
