// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_ASSOCIATIONS_H
#define LSST_JOINTCAL_ASSOCIATIONS_H

#include <string>
#include <iostream>
#include <list>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/afw/geom/Box.h"

#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/JointcalControl.h"

#include "lsst/afw/table/SortedCatalog.h"

namespace lsst {
namespace jointcal {

//! The class that implements the relations between MeasuredStar and FittedStar.
class Associations {
public:
    CcdImageList ccdImageList;      // the catalog handlers
    RefStarList refStarList;        // e.g. GAIA or SDSS reference stars
    FittedStarList fittedStarList;  // stars that are going to be fitted

    Point _commonTangentPoint;

public:
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
                                        std::numeric_limits<double>::quiet_NaN())) {}

    /// No moves or copies: jointcal only ever needs one Associations object.
    Associations(Associations const &) = delete;
    Associations(Associations &&) = delete;
    Associations &operator=(Associations const &) = delete;
    Associations &operator=(Associations &&) = delete;

    /**
     * @brief      Sets a shared tangent point for all ccdImages.
     *
     * @param      commonTangentPoint  The common tangent point of all input images (decimal degrees).
     */
    void setCommonTangentPoint(lsst::afw::geom::Point2D const &commonTangentPoint);

    //! can be used to project sidereal coordinates related to the image set on a plane.
    Point getCommonTangentPoint() const { return _commonTangentPoint; }

    /**
     * @brief      Create a ccdImage from an exposure catalog and metadata, and add it to the list.
     *
     * @param[in]  catalog    The extracted source catalog, selected for good astrometric sources.
     * @param[in]  wcs        The exposure's original wcs
     * @param[in]  visitInfo  The exposure's visitInfo object
     * @param[in]  bbox       The bounding box of the exposure
     * @param[in]  filter     The exposure's filter
     * @param[in]  photoCalib The exposure's photometric calibration
     * @param[in]  visit      The visit identifier
     * @param[in]  ccd        The ccd identifier
     * @param[in]  control    The JointcalControl object
     */
    void addImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &catalog,
                  std::shared_ptr<lsst::afw::image::TanWcs> wcs,
                  std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo, lsst::afw::geom::Box2I const &bbox,
                  std::string const &filter, std::shared_ptr<afw::image::PhotoCalib> photoCalib,
                  std::shared_ptr<afw::cameraGeom::Detector> detector, int visit, int ccd,
                  std::shared_ptr<lsst::jointcal::JointcalControl> control);

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
     * @param      refFluxMap     fluxes per filter of corresponding refCat objects (can be empty)
     * @param      refFluxErrMap  flux errors per filter of corresponding refCat objects (can be empty)
     */
    void collectRefStars(lsst::afw::table::SortedCatalogT<lsst::afw::table::SimpleRecord> &refCat,
                         afw::geom::Angle matchCut, std::string const &fluxField,
                         std::map<std::string, std::vector<double>> const &refFluxMap,
                         std::map<std::string, std::vector<double>> const &refFluxErrMap);

    //! Sends back the fitted stars coordinates on the sky FittedStarsList::inTangentPlaneCoordinates keeps
    //! track of that.
    void deprojectFittedStars();

//! Set the color field of FittedStar 's from a colored catalog.
/* If Color is "g-i", then the color is assigned from columns "g" and "i" of the colored catalog. */
#ifdef TODO
    void setFittedStarColors(std::string const &dicStarListName, std::string const &color,
                             double matchCutArcSec);
#endif

    /**
     * @brief      Apply quality cuts on potential FittedStars
     *
     * @param[in]  minMeasurements  The minimum number of measuredStars for a FittedStar to be included.
     */
    void selectFittedStars(int minMeasurements);

    CcdImageList const &getCcdImageList() const { return ccdImageList; }

    //! Number of different bands in the input image list. Not implemented so far
    unsigned getNFilters() const { return _filterMap.size(); }

    // Return the bounding box in (ra, dec) coordinates containing the whole catalog
    const lsst::afw::geom::Box2D getRaDecBBox();

    /**
     * @brief      return the number of CcdImages with non-empty catalogs to-be-fit.
     */
    int nCcdImagesValidForFit() const {
        return std::count_if(
                ccdImageList.begin(), ccdImageList.end(),
                [](std::shared_ptr<CcdImage> const &item) { return item->getCatalogForFit().size() > 0; });
    }

private:
    void associateRefStars(double matchCutInArcsec, const Gtransfo *gtransfo);

    void assignMags();

    // Map from filter name to index in each refStar's _refFlux/_refFluxErr vector.
    std::unordered_map<std::string, std::size_t> _filterMap;
};

}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_ASSOCIATIONS_H
