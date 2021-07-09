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

#ifndef LSST_JOINTCAL_ASSOCIATIONS_H
#define LSST_JOINTCAL_ASSOCIATIONS_H

#include <string>
#include <iostream>
#include <list>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/geom/Box.h"
#include "lsst/sphgeom/Circle.h"

#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/JointcalControl.h"

#include "lsst/afw/table/SortedCatalog.h"

namespace lsst {
namespace jointcal {

using RefFluxMapType = std::map<std::string, std::vector<double>>;

//! The class that implements the relations between MeasuredStar and FittedStar.
class Associations {
public:
    CcdImageList ccdImageList;      // the catalog handlers
    RefStarList refStarList;        // e.g. GAIA or SDSS reference stars
    FittedStarList fittedStarList;  // stars that are going to be fitted

    // These are strictly speaking not needed anymore (after DM-4043),
    // but keeping them seems cleaner then exposing the lists themselves.
    size_t refStarListSize() { return refStarList.size(); }
    size_t fittedStarListSize() { return fittedStarList.size(); }

    /**
     * Source selection is performed in python, so Associations' constructor
     * only initializes a couple of variables.
     */
    Associations()
            : _commonTangentPoint(Point(std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN())),
              _maxMeasuredStars(0) {}

    /**
     * Create an Associations object from a pre-built list of ccdImages.
     *
     * This is primarily useful for tests that build their own ccdImageList, but it could be used to help
     * parallelize the creation of the ccdImages.
     *
     * @param imageList A pre-built ccdImage list.
     * @param epoch The julian epoch year to which all proper motion corrections should be made.
     */
    Associations(CcdImageList const &imageList, double epoch = 0)
            : ccdImageList(imageList),
              _commonTangentPoint(Point(std::numeric_limits<double>::quiet_NaN(),
                                        std::numeric_limits<double>::quiet_NaN())),
              _maxMeasuredStars(0),
              _epoch(epoch) {}

    /// No moves or copies: jointcal only ever needs one Associations object.
    Associations(Associations const &) = delete;
    Associations(Associations &&) = delete;
    Associations &operator=(Associations const &) = delete;
    Associations &operator=(Associations &&) = delete;

    /**
     * Sets a shared tangent point for all ccdImages, using the mean of the centers of all ccdImages.
     */
    void computeCommonTangentPoint();

    /**
     * Shared tangent point for all ccdImages (decimal degrees).
     *
     * Used to project sidereal coordinates related to the images onto a single plane.
     */
    ///@{
    void setCommonTangentPoint(lsst::geom::Point2D const &commonTangentPoint);
    Point getCommonTangentPoint() const { return _commonTangentPoint; }
    ///@}

    /// The number of MeasuredStars at the start of fitting, before any outliers are removed.
    size_t getMaxMeasuredStars() const { return _maxMeasuredStars; }

    /**
     * Common epoch of all of the ccdImages as a Julian Epoch Year (e.g. 2000.0 for J2000).
     */
    ///@{
    double getEpoch() const { return _epoch; }
    void setEpoch(double epoch) { _epoch = epoch; }
    ///@}

    /**
     * @brief      Create a ccdImage from an exposure catalog and metadata, and add it to the list.
     *
     * @param[in]  catalog    The extracted source catalog, selected for good astrometric sources.
     * @param[in]  wcs        The exposure's original wcs
     * @param[in]  visitInfo  The exposure's visitInfo object
     * @param[in]  bbox       The bounding box of the exposure
     * @param[in]  filter     The exposure's filter
     * @param[in]  photoCalib The exposure's photometric calibration
     * @param[in]  detector   The exposure's detector
     * @param[in]  visit      The visit identifier
     * @param[in]  ccd        The ccd identifier
     * @param[in]  control    The JointcalControl object
     */
    void createCcdImage(afw::table::SourceCatalog &catalog, std::shared_ptr<lsst::afw::geom::SkyWcs> wcs,
                        std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, lsst::geom::Box2I const &bbox,
                        std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib,
                        std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
                        lsst::jointcal::JointcalControl const &control);

    /**
     * Add a pre-constructed ccdImage to the ccdImageList.
     */
    void addCcdImage(std::shared_ptr<CcdImage> const ccdImage) { ccdImageList.push_back(ccdImage); }

    //! incrementaly builds a merged catalog of all image catalogs
    void associateCatalogs(const double matchCutInArcsec = 0, const bool useFittedList = false,
                           const bool enlargeFittedList = true);

    /**
     * @brief      Collect stars from an external reference catalog and associate them with fittedStars.
     *
     * @param      refCat         The catalog of reference sources
     * @param[in]  matchCut       Separation radius to match fitted and
     *                            reference stars.
     * @param      fluxField      The field name in refCat to get the flux from.
     * @param      refCoordinateErr Error on reference catalog coordinates [mas]. If not NaN, this
     *                              overrides the `coord_*_err` values in the reference catalog itself.
     *                              This value is divided by cos(dec) before being used for ra_err.
     * @param      rejectBadFluxes  Reject reference sources with flux=NaN or 0 and/or fluxErr=NaN or 0.
     *                              Typically false for astrometry and true for photometry.
     */
    void collectRefStars(afw::table::SimpleCatalog &refCat, geom::Angle matchCut,
                         std::string const &fluxField, float refCoordinateErr, bool rejectBadFluxes = false);

    //! Sends back the fitted stars coordinates on the sky FittedStarsList::inTangentPlaneCoordinates keeps
    //! track of that.
    void deprojectFittedStars();

    /**
     * Prepare the fittedStar list by making quality cuts and normalizing measurements.
     *
     * @param[in]  minMeasurements  The minimum number of measuredStars for a FittedStar to be included.
     */
    void prepareFittedStars(int minMeasurements);

    /**
     * Remove FittedStars that have no measured stars; this can happen after outlier rejection.
     *
     * Use this to perform e.g. position minimization with outlier rejection after model minimization has
     * removed measuredStar outliers, to prevent the matrix from becoming non-positive definite.
     * Once this has been called, prepareFittedStars() has to be called again if the full
     * measuredStar->FittedStar relationship needs to be reconstructed.
     */
    void cleanFittedStars();

    CcdImageList const &getCcdImageList() const { return ccdImageList; }

    //! Number of different bands in the input image list. Not implemented so far
    unsigned getNFilters() const { return 1; }

    /**
     * Return the bounding circle in on-sky (RA, Dec) coordinates containing all CcdImages.
     *
     * Requires that computeCommonTangentPoint() be called first, so that sensor bounding boxes can be
     * transformed into the common tangent plane.
     */
    lsst::sphgeom::Circle computeBoundingCircle() const;

    /**
     * @brief      return the number of CcdImages with non-empty catalogs to-be-fit.
     */
    int nCcdImagesValidForFit() const;

    /**
     * @brief      Return the number of fittedStars that have an associated refStar.
     */
    size_t nFittedStarsWithAssociatedRefStar() const;

private:
    void associateRefStars(double matchCutInArcsec, const AstrometryTransform *transform);

    void assignMags();

    /**
     * Apply quality cuts on potential FittedStars
     *
     * @param[in]  minMeasurements  The minimum number of measuredStars for a FittedStar to be included.
     */
    void selectFittedStars(int minMeasurements);

    /**
     * Make fitted star positions and fluxes be the average of their measured stars.
     *
     * Only call after selectFittedStars() has been called: it assumes that each measuredStar points to a
     * fittedStar, and that the measurementCount for each fittedStar is correct.
     */
    void normalizeFittedStars();

    // Common tangent point on-sky of all of the ccdImages, typically determined by computeCommonTangentPoint.
    // (decimal degrees)
    Point _commonTangentPoint;

    // The number of MeasuredStars at the start of fitting, before any outliers are removed.
    // This is used to reserve space in vectors for e.g. outlier removal, but is not updated during outlier
    // removal or cleanup, so should only be used as an upper bound on the number of MeasuredStars.
    size_t _maxMeasuredStars;

    // Julian Epoch Year (e.g. 2000.0 for J2000)
    // Common epoch of all of the ccdImages, typically computed externally via astropy and then set.
    double _epoch;
};

}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_ASSOCIATIONS_H
