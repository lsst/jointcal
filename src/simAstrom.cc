#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <stdlib.h> /* for getenv */

#define _GNU_SOURCE 1
#define __USE_GNU
#include <fenv.h> 

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
#include "lsst/meas/simastrom/simAstrom.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Point.h"
#include "lsst/meas/simastrom/Associations.h"
#include "lsst/meas/simastrom/Projectionhandler.h"
#include "lsst/meas/simastrom/SimplePolyModel.h"
#include "lsst/meas/simastrom/AstromFit.h"

#include "Eigen/Core"

namespace pexExcept = lsst::pex::exceptions;
namespace afwTable = lsst::afw::table;

namespace lsst {
namespace meas {
namespace simastrom {
    
    static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked. */

  if (getenv("DUMP_CORE_ON_FPE"))
    feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
//    feenableexcept (FE_ALL_EXCEPT);
} 
    
    simAstrom::simAstrom(
        std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> > const sourceList,
        std::vector<PTR(lsst::daf::base::PropertySet)> const metaList,
        std::vector<PTR(lsst::afw::image::TanWcs)> const wcsList,
        std::vector<lsst::afw::geom::Box2I> const bboxList,
        std::vector<std::string> const filterList,
        std::vector<PTR(lsst::afw::image::Calib)> const calibList,
        std::vector<int> const visitList,
        std::vector<int> const ccdList
    ):
        _sourceList(sourceList),
        _metaList(metaList),
        _wcsList(wcsList),
        _bboxList(bboxList),
        _filterList(filterList),
        _calibList(calibList),
        _visitList(visitList),
        _ccdList(ccdList)
{
    std::cout << "simAstrom constructor invoked " << std::endl;
    std::cout << "Vectors contain : " << _sourceList.size() << " elements" << std::endl;
    std::cout << _metaList[0]->get<double>("LATITUDE") << std::endl;
    std::cout << _sourceList[1][10].getRa() << std::endl;
    std::cout << _wcsList[1]->getPixelOrigin() << std::endl;
    std::cout << _bboxList[1] << std::endl;
    std::cout << _filterList[1] << std::endl;
    std::cout << _visitList[1] << std::endl;
    std::cout << _ccdList[1] << std::endl;
    
    // Check how to get SIP coefficients from WCS
    lsst::daf::base::PropertyList::Ptr wcsMeta = _wcsList[1]->getFitsMetadata();
//    std::cout << wcsMeta->getOrderedNames() << std::endl;
    std::cout << wcsMeta->get<int>("A_ORDER") << std::endl;
    
    Eigen::MatrixXd sipA;
    lsst:afw::image::TanWcs::decodeSipHeader(*wcsMeta, "A", sipA);
    std::cout << sipA << std::endl;
    
    //    std::cout << _sourceList[1].getSchema() << std::endl;
//    auto centroidKey = _sourceList[1].getSchema().find<double>("base_SdssCentroid_x").key;
//    for (afwTable::SourceCatalog::const_iterator sourcePtr = _sourceList[1].begin();
//            sourcePtr != _sourceList[1].end(); ++sourcePtr) {
//                std::cout << sourcePtr->getCentroid() << std::endl;
//                std::cout << sourcePtr->get(centroidKey) << std::endl;
//            } 
    
    // Create and load an Associations object
    Associations *assoc = new Associations();    
    for (int i=0; i<_sourceList.size(); i++) {
        assoc->AddImage(_sourceList[i], _wcsList[i], _metaList[i], _bboxList[i], _filterList[i], 
        _calibList[i], _visitList[i], _ccdList[i]);
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

    AstromFit astromFit(*assoc, &spm); 

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

      for (unsigned k=0; k<5; ++k)
        {
          astromFit.RemoveOutliers(5.);
          std::cout << "After outliers removal" << std::endl;
          astromFit.Minimize("Positions Distortions");
          std::cout << astromFit.ComputeChi2() << std::endl;
        }
      astromFit.MakeResTuple("res.list");
}    
    
}}}
