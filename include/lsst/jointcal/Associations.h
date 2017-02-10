// -*- C++ -*-
//
#ifndef ASSOCIATIONS__H
#define ASSOCIATIONS__H

#include <string>
#include <iostream>

//#include "stringlist.h"

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
#include "lsst/jointcal/Jointcal.h"

#include "lsst/afw/table/SortedCatalog.h"

namespace lsst {
namespace jointcal {

//! The class that implements the relations between MeasuredStar and FittedStar.
class Associations {
public:

    CcdImageList ccdImageList; // the catalog handlers
    RefStarList refStarList; // e.g. GAIA or SDSS reference stars
    FittedStarList fittedStarList; // stars that are going to be fitted

    // TODO: NOTE: these only exist because those lists can't be swig'd because
    // swig doesn't like boost::intrusive_ptr (which the lists contain).
    // Once DM-4043 is solved, delete these.
    size_t refStarListSize()
    { return refStarList.size(); }
    size_t fittedStarListSize()
    { return fittedStarList.size(); }

    // fit cuts and stuff:
    Point _commonTangentPoint;

    void assignMags();

public:

    /**
     * Source selection is performed in python, so Associations' constructor
     * only initializes a couple of variables.
     */
    Associations();

    //! Sets a tangent point (reasonably centered for the input image set).
    void setCommonTangentPoint(Point const &commonTangentPoint)
    { _commonTangentPoint = commonTangentPoint;};

    //! can be used to project sidereal coordinates related to the image set on a plane.
    Point getCommonTangentPoint() const { return _commonTangentPoint;}

    //! same as above
    bool addImage(lsst::afw::table::SortedCatalogT<lsst::afw::table::SourceRecord> &catalog,
                  std::shared_ptr<lsst::afw::image::TanWcs> wcs,
                  std::shared_ptr<lsst::afw::image::VisitInfo> visitInfo,
                  lsst::afw::geom::Box2I const &bbox,
                  std::string const &filter,
                  std::shared_ptr<lsst::afw::image::Calib> calib,
                  int visit,
                  int ccd,
                  std::shared_ptr<lsst::jointcal::JointcalControl> control);

    //! incrementaly builds a merged catalog of all image catalogs
    void associateCatalogs(const double matchCutInArcsec = 0,
                           const bool UseFittedList = false,
                           const bool EnlargeFittedList = true);

    //! Collect stars from an external reference catalog
    void collectRefStars(lsst::afw::table::SortedCatalogT< lsst::afw::table::SimpleRecord > &Ref,
                         std::string const &fluxField);

    //! Sends back the fitted stars coordinates on the sky FittedStarsList::inTangentPlaneCoordinates keeps track of that.
    void deprojectFittedStars();


    //! Set the color field of FittedStar 's from a colored catalog.
    /* If Color is "g-i", then the color is assigned from columns "g" and "i" of the colored catalog. */
#ifdef TODO
    void setFittedStarColors(std::string const &DicStarListName,
                             std::string const &Color,
                             double MatchCutArcSec);
#endif

    //! apply cuts (mainly number of measurements) on potential FittedStars
    void selectFittedStars();

    const CcdImageList& getCcdImageList() const {return ccdImageList;}

    //! Number of different bands in the input image list. Not implemented so far
    unsigned NBands() const {return 1;}

    // Return the bounding box in (ra, dec) coordinates containing the whole catalog
    const lsst::afw::geom::Box2D getRaDecBBox();


private:
    void associateRefStars(double matchCutInArcsec, const Gtransfo* gtransfo);

};

#endif /* ASSOCIATIONS__H */

}
} // end of namespaces
