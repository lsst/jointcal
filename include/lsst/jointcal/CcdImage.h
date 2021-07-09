// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_JOINTCAL_CCD_IMAGE_H
#define LSST_JOINTCAL_CCD_IMAGE_H

#include <list>
#include <string>

#include "lsst/afw/cameraGeom/Detector.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/geom/Box.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Frame.h"

namespace lsst {
namespace jointcal {

typedef std::list<std::shared_ptr<CcdImage>> CcdImageList;

typedef int VisitIdType;
typedef int CcdIdType;
/// For hashing a ccdImage: the pair of (visit, ccd) IDs should be unique to each ccdImage.
struct CcdImageKey {
    VisitIdType visit;
    CcdIdType ccd;
    bool operator!=(CcdImageKey const &right) const { return !(*this == right); }
    bool operator==(CcdImageKey const &right) const { return (visit == right.visit) && (ccd == right.ccd); }
};
// typedef std::pair<VisitIdType, CcdIdType> CcdImageKey;
std::ostream &operator<<(std::ostream &out, CcdImageKey const &key);

/**
 * Handler of an actual image from a single CCD.
 * NOTE: could possibly be replaced with a subclass of afw.image.Exposure?
 */
class CcdImage {
public:
    CcdImage(afw::table::SourceCatalog &record, std::shared_ptr<lsst::afw::geom::SkyWcs> wcs,
             std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, lsst::geom::Box2I const &bbox,
             std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib,
             std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
             std::string const &fluxField);

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

    /// Clear the catalog for fitting and set it to a copy of the whole catalog.
    void resetCatalogForFit() {
        getCatalogForFit().clear();
        getWholeCatalog().copyTo(getCatalogForFit());
    }

    /**
     * Count the number of valid measured and reference stars that fall within this ccdImage.
     *
     * Measured stars are counted if they are valid. Reference stars are counted if a measured star
     * (valid or not) has a fittedStar that has an associated refStar.
     *
     * @return Number of (measured, reference) stars in the image.
     */
    std::pair<int, int> countStars() const;

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

    std::shared_ptr<AstrometryTransform> const getPixelToCommonTangentPlane() const {
        return _pixelToCommonTangentPlane;
    }

    std::shared_ptr<AstrometryTransform> const getCommonTangentPlaneToTangentPlane() const {
        return _commonTangentPlaneToTangentPlane;
    }

    std::shared_ptr<AstrometryTransform> const getTangentPlaneToCommonTangentPlane() const {
        return _tangentPlaneToCommonTangentPlane;
    }

    std::shared_ptr<AstrometryTransform> const getPixelToTangentPlane() const { return _pixelToTangentPlane; }

    std::shared_ptr<AstrometryTransform> const getSkyToTangentPlane() const { return _skyToTangentPlane; }

    //! returns ccd ID
    CcdIdType getCcdId() const { return _ccdId; }

    //! returns visit ID
    VisitIdType getVisit() const { return _visit; }

    std::shared_ptr<afw::cameraGeom::Detector> getDetector() const { return _detector; }

    CcdImageKey getHashKey() const { return CcdImageKey{_visit, _ccdId}; }

    //!  Airmass
    double getAirMass() const { return _airMass; }

    double getEpoch() const { return _epoch; }

    //! Return the exposure's photometric calibration
    std::shared_ptr<afw::image::PhotoCalib> getPhotoCalib() const { return _photoCalib; }

    /**
     * @brief      Gets the boresight RA/Dec.
     */
    lsst::geom::SpherePoint getBoresightRaDec() const { return _boresightRaDec; }

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
    std::shared_ptr<AstrometryTransform> const getReadWcs() const { return _readWcs; }

    //! Frame in pixels
    Frame const &getImageFrame() const { return _imageFrame; }

private:
    void loadCatalog(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> const &Cat,
                     std::string const &fluxField);

    jointcal::Frame _imageFrame;  // in pixels

    MeasuredStarList _wholeCatalog;  // the catalog of measured objets
    MeasuredStarList _catalogForFit;

    std::shared_ptr<AstrometryTransformSkyWcs> _readWcs;  // apply goes from pix to sky

    // The following ones should probably be mostly removed.
    // go from CommonTangentPlane to this tangent plane.
    std::shared_ptr<AstrometryTransform> _commonTangentPlaneToTangentPlane;
    std::shared_ptr<AstrometryTransform> _tangentPlaneToCommonTangentPlane;  // reverse one
    std::shared_ptr<AstrometryTransform> _pixelToCommonTangentPlane;         // pixels -> CTP
    std::shared_ptr<AstrometryTransform> _pixelToTangentPlane;

    std::shared_ptr<AstrometryTransform> _skyToTangentPlane;

    std::string _name;
    CcdIdType _ccdId;
    VisitIdType _visit;

    lsst::geom::SpherePoint _boresightRaDec;
    double _airMass;  // airmass value.
    double _epoch;    // julian epoch year (e.g. 2000.0 for J2000)
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
    size_t operator()(lsst::jointcal::CcdImageKey const &key) const {
        return hash<size_t>()(static_cast<size_t>(key.visit) | (static_cast<size_t>(key.ccd) << 32));
    }
};
}  // namespace std

#endif  // LSST_JOINTCAL_CCD_IMAGE_H
