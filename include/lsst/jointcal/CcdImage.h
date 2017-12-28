// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CCD_IMAGE_H
#define LSST_JOINTCAL_CCD_IMAGE_H

#include <list>
#include <string>

#include "lsst/afw/cameraGeom/Detector.h"
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
/// For hashing a ccdImage: the pair of (visit, ccd) IDs should be unique to each ccdImage.
typedef std::pair<VisitIdType, CcdIdType> CcdImageKey;
std::ostream &operator<<(std::ostream &out, CcdImageKey const &key);

/**
 * Handler of an actual image from a single CCD.
 * NOTE: could possibly be replaced with a subclass of afw.image.Exposure?
 */
class CcdImage {
public:
    CcdImage(afw::table::SourceCatalog &record, std::shared_ptr<lsst::afw::image::TanWcs> wcs,
             std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, afw::geom::Box2I const &bbox,
             std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib,
             std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
             std::string const &fluxField);

    // For pickling
    CcdImage(jointcal::Frame imageFrame, MeasuredStarList wholeCatalog, MeasuredStarList catalogForFit,
             std::shared_ptr<BaseTanWcs> readWcs, std::shared_ptr<Gtransfo> inverseReadWcs,
             std::shared_ptr<Gtransfo> commonTangentPlane2TP, std::shared_ptr<Gtransfo> TP2CommonTangentPlane,
             std::shared_ptr<Gtransfo> pix2CommonTangentPlane, std::shared_ptr<Gtransfo> pix2TangentPlane,
             std::shared_ptr<Gtransfo> sky2TP, std::string name, CcdIdType ccdId, VisitIdType visit,
             afw::coord::IcrsCoord boresightRaDec, double airMass, double mjd,
             std::shared_ptr<afw::image::PhotoCalib> photoCalib,
             std::shared_ptr<afw::cameraGeom::Detector> detector, double sinEta, double cosEta, double tanZ,
             double lstObs, double hourAngle, std::string filter, jointcal::Point commonTangentPoint)
            : _imageFrame(imageFrame),
              _wholeCatalog(wholeCatalog),
              _catalogForFit(catalogForFit),
              _readWcs(readWcs),
              _inverseReadWcs(inverseReadWcs),
              _commonTangentPlane2TP(commonTangentPlane2TP),
              _TP2CommonTangentPlane(TP2CommonTangentPlane),
              _pix2CommonTangentPlane(pix2CommonTangentPlane),
              _pix2TangentPlane(pix2TangentPlane),
              _sky2TP(sky2TP),
              _name(name),
              _ccdId(ccdId),
              _visit(visit),
              _boresightRaDec(boresightRaDec),
              _airMass(airMass),
              _mjd(mjd),
              _photoCalib(photoCalib),
              _detector(detector),
              _sinEta(sinEta),
              _cosEta(cosEta),
              _tanZ(tanZ),
              _lstObs(lstObs),
              _hourAngle(hourAngle),
              _filter(filter),
              _commonTangentPoint(commonTangentPoint) {}

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
    jointcal::Point const &getCommonTangentPoint() const { return _commonTangentPoint; }

    std::shared_ptr<Gtransfo> const getPix2CommonTangentPlane() const { return _pix2CommonTangentPlane; }

    std::shared_ptr<Gtransfo> const getCommonTangentPlane2TP() const { return _commonTangentPlane2TP; }

    std::shared_ptr<Gtransfo> const getTP2CommonTangentPlane() const { return _TP2CommonTangentPlane; }

    std::shared_ptr<Gtransfo> const getPix2TangentPlane() const { return _pix2TangentPlane; }

    std::shared_ptr<Gtransfo> const getSky2TP() const { return _sky2TP; }

    //! returns ccd ID
    CcdIdType getCcdId() const { return _ccdId; }

    //! returns visit ID
    VisitIdType getVisit() const { return _visit; }

    std::shared_ptr<afw::cameraGeom::Detector> getDetector() const { return _detector; }

    CcdImageKey getHashKey() const { return CcdImageKey(_visit, _ccdId); }

    //!  Airmass
    double getAirMass() const { return _airMass; }

    //! Julian Date
    double getMjd() const { return _mjd; }

    //! Return the exposure's photometric calibration
    std::shared_ptr<afw::image::PhotoCalib> getPhotoCalib() const { return _photoCalib; }

    /**
     * @brief      Gets the boresight RA/Dec.
     */
    lsst::afw::coord::IcrsCoord getBoresightRaDec() const { return _boresightRaDec; }

    double getHourAngle() const { return _hourAngle; }

    double getLstObs() const { return _lstObs; }

    //! Parallactic angle
    double getSinEta() const { return _sinEta; }

    //! Parallactic angle
    double getCosEta() const { return _cosEta; }

    //! Parallactic angle
    double getTanZ() const { return _tanZ; }

    //!
    Point getRefractionVector() const { return Point(_tanZ * _cosEta, _tanZ * _sinEta); }

    //! return the CcdImage filter name
    std::string getFilter() const { return _filter; }

    //! the wcs read in the header. NOT updated when fitting.
    std::shared_ptr<Gtransfo> const getReadWcs() const { return _readWcs; }

    //! the inverse of the one above.
    std::shared_ptr<Gtransfo> const getInverseReadWcs() const { return _inverseReadWcs; }

    //! Frame in pixels
    Frame const &getImageFrame() const { return _imageFrame; }

private:
    void loadCatalog(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> const &Cat,
                     std::string const &fluxField);

    jointcal::Frame _imageFrame;  // in pixels

    MeasuredStarList _wholeCatalog;  // the catalog of measured objets
    MeasuredStarList _catalogForFit;

    std::shared_ptr<BaseTanWcs> _readWcs;       // i.e. from pix to sky
    std::shared_ptr<Gtransfo> _inverseReadWcs;  // i.e. from sky to pix

    // The following ones should probably be mostly removed.
    // go from CommonTangentPlane to this tangent plane.
    std::shared_ptr<Gtransfo> _commonTangentPlane2TP;
    std::shared_ptr<Gtransfo> _TP2CommonTangentPlane;   // reverse one
    std::shared_ptr<Gtransfo> _pix2CommonTangentPlane;  // pixels -> CTP
    std::shared_ptr<Gtransfo> _pix2TangentPlane;

    std::shared_ptr<Gtransfo> _sky2TP;

    std::string _name;
    CcdIdType _ccdId;
    VisitIdType _visit;

    lsst::afw::coord::IcrsCoord _boresightRaDec;
    double _airMass;  // airmass value.
    double _mjd;      // modified julian date
    std::shared_ptr<afw::image::PhotoCalib> _photoCalib;
    std::shared_ptr<afw::cameraGeom::Detector> _detector;
    // refraction
    // eta : parallactic angle, z: zenithal angle (X = 1/cos(z))
    double _sinEta, _cosEta, _tanZ;
    // Local Sidereal Time and hour angle of observation
    double _lstObs, _hourAngle;

    std::string _filter;

    jointcal::Point _commonTangentPoint;
};
}  // namespace jointcal
}  // namespace lsst

// Add our preferred hash of CcdImageKey to the std:: namespace, so it's always available "for free".
namespace std {
template <>
/**
 * Hash a ccdImage by its visit and ccd IDs.
 *
 * ccdId and visitId are both 32-bit ints, hash() returns a size_t, so put the ccdId in the
 * most-significant-bit, and the visitId in the least for a simple, unique, hash per ccdImage.
 */
struct hash<lsst::jointcal::CcdImageKey> {
    size_t operator()(lsst::jointcal::CcdImageKey const &ccdImage) const {
        return hash<size_t>()(static_cast<size_t>(ccdImage.first) |
                              (static_cast<size_t>(ccdImage.second) << 32));
    }
};
}  // namespace std

#endif  // LSST_JOINTCAL_CCD_IMAGE_H
