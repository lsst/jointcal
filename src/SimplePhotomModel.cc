#include <iostream>

#include "lsst/jointcal/SimplePhotomModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace lsst {
namespace jointcal {

SimplePhotomModel::SimplePhotomModel(const CcdImageList &ccdImageList)
{
  unsigned refVisit = -1;
  for (auto i = ccdImageList.begin(); i !=ccdImageList.end(); ++i)
    {
      const CcdImage &im = **i;
      unsigned visit = im.getVisit();
      if (refVisit == -1) refVisit = visit;
      if (visit==refVisit)
	  _myMap[&im].fixed=true;
      else _myMap[&im].fixed=false;
    }
  std::cout << "INFO: SimplePhotomModel : using exposure " << refVisit << " as photometric reference " << std::endl;
}

unsigned SimplePhotomModel::assignIndices(const std::string &whatToFit,  unsigned firstIndex)
{
  unsigned ipar = firstIndex;
  for (auto i = _myMap.begin(); i!= _myMap.end(); ++i)
    {
      PhotomStuff& pf=i->second;
      if (pf.fixed) continue;
      pf.index=ipar;
      ipar++;
    }
  return ipar;
}

 void SimplePhotomModel::offsetParams(const Eigen::VectorXd &delta)
 {
   for (auto i = _myMap.begin(); i!= _myMap.end(); ++i)
     {
       PhotomStuff& pf=i->second;
       if (!pf.fixed) pf.factor += delta[pf.index];
     }
 }

 SimplePhotomModel::PhotomStuff& SimplePhotomModel::find(const CcdImage &ccdImage)
   {
     auto i = _myMap.find(&ccdImage);
     if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::find, cannot find CcdImage "+ccdImage.getName());
     return (i->second);
   }

 const SimplePhotomModel::PhotomStuff& SimplePhotomModel::find(const CcdImage &ccdImage)  const
   {
     auto i = _myMap.find(&ccdImage);
     if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::find, cannot find CcdImage "+ccdImage.getName());
     return (i->second);
   }


  double SimplePhotomModel::photomFactor(const CcdImage &ccdImage, const Point &where) const
 {
   const PhotomStuff &pf = find(ccdImage);
   return pf.factor;
 }

 void SimplePhotomModel::getIndicesAndDerivatives(const MeasuredStar &measuredStar,
                                                  const CcdImage &ccdImage,
                                                  std::vector<unsigned> &indices,
                                                  Eigen::VectorXd &D)
 {
   PhotomStuff &pf = find(ccdImage);
   if (pf.fixed) {indices.resize(0); return;}
   indices.resize(1);
   indices[0] = pf.index;
   D[0] = 1;
 }

}} // end of namespaces
