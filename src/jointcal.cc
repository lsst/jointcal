#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <stdlib.h> /* for getenv */

//#include "boost/scoped_array.hpp"
//#include "boost/shared_array.hpp"
//#include "boost/multi_index_container.hpp"
//#include "boost/multi_index/sequenced_index.hpp"
//#include "boost/multi_index/ordered_index.hpp"
//#include "boost/multi_index/global_fun.hpp"

#include "boost/make_shared.hpp"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/jointcal/jointcal.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/AstromFit.h"

#include "Eigen/Core"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace jointcal {
    
    void JointcalControl::validate() const {
        if (sourceFluxField.empty()) {
            throw LSST_EXCEPT(pexExcept::InvalidParameterError, "sourceFluxField must be specified");
        }
        std::cout << sourceFluxField << std::endl;
    }
    
    Jointcal::Jointcal(
        std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > const sourceList,
        std::vector<PTR(lsst::daf::base::PropertySet)> const metaList,
        std::vector<PTR(lsst::afw::image::TanWcs)> const wcsList,
        std::vector<lsst::afw::geom::Box2I> const bboxList,
        std::vector<std::string> const filterList,
        std::vector<PTR(lsst::afw::image::Calib)> const calibList,
        std::vector<int> const visitList,
        std::vector<int> const ccdList,
        std::vector<std::string> const cameraList,
        PTR(lsst::jointcal::JointcalControl) const control
    ):
        _sourceList(sourceList),
        _metaList(metaList),
        _wcsList(wcsList),
        _bboxList(bboxList),
        _filterList(filterList),
        _calibList(calibList),
        _visitList(visitList),
        _ccdList(ccdList),
        _cameraList(cameraList)
{

    std::cout << "sourceFluxField is set to : " << control->sourceFluxField << std::endl;
    if (_wcsList.size() == 0)
      {
	std::cout << "jointcal::Jointcal : empty image list, we give up" << std::endl;
	return;
      }

    
    // Check how to get SIP coefficients from WCS
    lsst::daf::base::PropertyList::Ptr wcsMeta = _wcsList[0]->getFitsMetadata();
//    std::cout << wcsMeta->getOrderedNames() << std::endl;
    std::cout << wcsMeta->get<int>("A_ORDER") << std::endl;
    
    Eigen::MatrixXd sipA;
    afw::image::TanWcs::decodeSipHeader(*wcsMeta, "A", sipA);
    std::cout << sipA << std::endl;
    
    //    std::cout << _sourceList[1].getSchema() << std::endl;
//    auto centroidKey = _sourceList[1].getSchema().find<double>("base_SdssCentroid_x").key;
//    for (afwTable::SourceCatalog::const_iterator sourcePtr = _sourceList[1].begin();
//            sourcePtr != _sourceList[1].end(); ++sourcePtr) {
//                std::cout << sourcePtr->getCentroid() << std::endl;
//                std::cout << sourcePtr->get(centroidKey) << std::endl;
//            }
    
    // Create and load an Association object
    Associations *assoc = new Associations();
    for (unsigned i=0; i<_sourceList.size(); i++) {
        assoc->AddImage(_sourceList[i], _wcsList[i], _metaList[i], _bboxList[i], _filterList[i],
        _calibList[i], _visitList[i], _ccdList[i], _cameraList[i], control);
    }
    
    // Associates catalog
    double matchCut = 3.0;   // Should be passed as a parameter
    assoc->AssociateCatalogs(matchCut);
    assoc->CollectRefStars(/*ProjectOnTP =  */ false);
    assoc->SelectFittedStars();
    
    TanPix2RaDec ctp2Sky(GtransfoLin(), assoc->CommonTangentPoint());
    
    assoc->fittedStarList.ApplyTransfo(ctp2Sky);
    OneTPPerShoot sky2TP(assoc->TheCcdImageList());
    SimplePolyModel spm(assoc->TheCcdImageList(), &sky2TP, true, 0);

    double posError = 0.02;
    AstromFit astromFit(*assoc, &spm, posError);

    std::string whatToFit = "Distortions";

    for (unsigned iter=0; iter<2;++iter)
      {
	std::cout << " Fitting only mappings" << std::endl;
	std::cout << astromFit.ComputeChi2() << std::endl;
	astromFit.Minimize(whatToFit);
      }
      
    whatToFit = "Positions";
    for (unsigned iter=0; iter<2;++iter)
        {
        std::cout << " Fitting only positions" << std::endl;
        std::cout << astromFit.ComputeChi2() << std::endl;
        astromFit.Minimize(whatToFit);
        }

    std::cout << astromFit.ComputeChi2() << std::endl;
    
    std::cout << " fitting positions and mappings" << std::endl;
      
    astromFit.Minimize("Positions Distortions");

    std::cout << astromFit.ComputeChi2() << std::endl;

    astromFit.MakeResTuple("res0.list");

    for (unsigned k=0; k<20; ++k)
        {
          astromFit.RemoveOutliers(5.,"Meas Ref");
          std::cout << "After outliers removal" << std::endl;
          astromFit.Minimize("Positions Distortions");
          std::cout << astromFit.ComputeChi2() << std::endl;
        }
    astromFit.MakeResTuple("res.list");
}
    
}} // end of namespaces
