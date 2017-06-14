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

typedef std::list<std::shared_ptr<CcdImage> > CcdImageList;

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

    void LoadCatalog(const lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &Cat,
                     const std::string &fluxField);

public:
    CcdImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &record,
             const PTR(lsst::afw::image::TanWcs) wcs, const PTR(lsst::afw::image::VisitInfo) visitInfo,
             const lsst::afw::geom::Box2I &bbox, const std::string &filter,
             const std::shared_ptr<afw::image::PhotoCalib> photoCalib, const int &visit, const int &ccd,
             const std::string &fluxField);

    //! Return the _name that identifies this ccdImage.
    std::string getName() const { return _name; }

    /**
     * @brief      Gets the as-read catalog.
     *
     * @return     The whole catalog.
     */
    const MeasuredStarList &getWholeCatalog() const { return _wholeCatalog; }

    //@{
    /**
     * @brief      Gets the catalog to be used for fitting, which may have been cleaned-up.
     *
     * @return     The catalog for fitting.
     */
    const MeasuredStarList &getCatalogForFit() const { return _catalogForFit; }
    MeasuredStarList &getCatalogForFit() { return _catalogForFit; }
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
    Point const &getCommonTangentPoint() const { return _commonTangentPoint; }

    //!
    const Gtransfo *getPix2CommonTangentPlane() const { return _pix2CommonTangentPlane.get(); }

    //!
    const Gtransfo *getCommonTangentPlane2TP() const { return _CTP2TP.get(); }

    //!
    const Gtransfo *getTP2CommonTangentPlane() const { return _TP2CTP.get(); }

    //!
    const Gtransfo *getPix2TangentPlane() const { return _pix2TP.get(); }

    //!
    const Gtransfo *getSky2TP() const { return _sky2TP.get(); }

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
    const Gtransfo *readWCS() const { return _readWcs.get(); }

    //! the inverse of the one above.
    const Gtransfo *getInverseReadWCS() const { return _inverseReadWcs.get(); }

    //! Frame in pixels
    const Frame &getImageFrame() const { return _imageFrame; }

private:
    CcdImage(const CcdImage &);  // forbid copies
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_CCD_IMAGE_H
