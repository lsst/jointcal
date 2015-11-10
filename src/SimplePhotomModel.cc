#include <iostream>

#include "lsst/meas/simastrom/SimplePhotomModel.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/MeasuredStar.h"

namespace lsst {
namespace meas {
namespace simastrom {



SimplePhotomModel::SimplePhotomModel(const CcdImageList &L)
{
  unsigned refShoot = -1;
  for (auto i = L.begin(); i !=L.end(); ++i)
    {
      const CcdImage &im = **i;
      unsigned shoot = im.Shoot();
      if (refShoot == -1) refShoot = shoot;
      if (shoot==refShoot)
	  _myMap[&im].fixed=true;
      else _myMap[&im].fixed=false;
    }
  std::cout << "INFO: SimplePhotomModel : using exposure " << refShoot << " as photometric reference " << std::endl;
}

unsigned SimplePhotomModel::AssignIndices(const std::string &WhatToFit,  unsigned FirstIndex)
{
  unsigned ipar = FirstIndex;
  for (auto i = _myMap.begin(); i!= _myMap.end(); ++i)
    {
      PhotomStuff& pf=i->second;
      if (pf.fixed) continue;
      pf.index=ipar;
      ipar++;
    }
  return ipar;
}

 void SimplePhotomModel::OffsetParams(const Eigen::VectorXd &Delta)
 {
   for (auto i = _myMap.begin(); i!= _myMap.end(); ++i)
     {
       PhotomStuff& pf=i->second;
       if (!pf.fixed) pf.factor += Delta[pf.index];
     }
 }

 SimplePhotomModel::PhotomStuff& SimplePhotomModel::find(const CcdImage &C) 
   {
     auto i = _myMap.find(&C);
     if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::find, never heard of CcdImage "+C.Name());
     return (i->second);
   }

 const SimplePhotomModel::PhotomStuff& SimplePhotomModel::find(const CcdImage &C)  const
   {
     auto i = _myMap.find(&C);
     if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::find, never heard of CcdImage "+C.Name());
     return (i->second);
   }


  double SimplePhotomModel::PhotomFactor(const CcdImage &C, const Point &Where) const
 {
   const PhotomStuff &pf = find(C);
   return pf.factor;
 }
     
 void SimplePhotomModel::GetIndicesAndDerivatives(const MeasuredStar &M,
						  const CcdImage &Ccd, 
						  std::vector<unsigned> &Indices,
						  Eigen::VectorXd &D)
 {
   PhotomStuff &pf = find(Ccd);
   if (pf.fixed) {Indices.resize(0); return;}
   Indices.resize(1);
   Indices[0] = pf.index;
   D[0] = 1;
 }



}}}
