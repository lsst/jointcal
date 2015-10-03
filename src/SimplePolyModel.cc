#include "lsst/meas/simastrom/Eigenstuff.h"
#include "lsst/meas/simastrom/SimplePolyModel.h"
#include "lsst/meas/simastrom/SimplePolyMapping.h"
#include "lsst/meas/simastrom/CcdImage.h"
#include "lsst/meas/simastrom/Projectionhandler.h"
#include "lsst/pex/exceptions.h"
#include <string>

#include "lsst/meas/simastrom/Gtransfo.h"

//const int distortionDegree=3;

namespace lsst {
namespace meas {
namespace simastrom {

// need a way to propagate the requested degree !
SimplePolyModel::SimplePolyModel(const CcdImageList &L, 
				 const ProjectionHandler* ProjH, 
				 bool InitFromWCS,
				 unsigned NNotFit,
				 unsigned degree) : _sky2TP(ProjH)

{
  // from datacards (or default)
//  unsigned degree = distortionDegree; 
  unsigned count = 0;
  for (auto i=L.cbegin(); i!= L.end(); ++i, ++count)
    {
      const CcdImage &im = **i;
      if (count < NNotFit)
	{
	  SimpleGtransfoMapping * id = new SimpleGtransfoMapping(GtransfoIdentity());
	  _myMap[&im] = id;
	  id->SetIndex(-1); // non sense, because it has no parameters
	}
      else
	// Given how AssignIndices works, only the SimplePolyMapping's
	// will actually be fitted, as NNotFit requests.
	{
	  GtransfoPoly pol(degree);
	  const Frame &frame  = im.ImageFrame();
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);
	  if (InitFromWCS) 
	    {	
	      pol = GtransfoPoly(im.Pix2TangentPlane(),
				 frame,
				 degree);
	      pol = pol*shiftAndNormalize.invert();

	    }
	  _myMap[&im] = new SimplePolyMapping(shiftAndNormalize, pol);
	}
    }
}


const Mapping* SimplePolyModel::GetMapping(const CcdImage &C) const
{
  mapType::const_iterator i = _myMap.find(&C);
  if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::GetMapping, never heard of CcdImage "+C.Name());
  return (i->second);
}

unsigned SimplePolyModel::AssignIndices(unsigned FirstIndex, std::string &WhatToFit)
{
  if (WhatToFit.find("Distortions") == std::string::npos)
    { 
      std::cout << "SimplePolyModel::AssignIndices is called and Distortions is *not*  in WhatToFit" << std::endl;
      return 0;
    }
  unsigned index = FirstIndex;
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    {
      SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
      if (!p) continue; // it should be GtransfoIdentity
      p->SetIndex(index);
      index+= p->Npar();
    }
  return index;
}

void SimplePolyModel::OffsetParams(const Eigen::VectorXd &Delta)
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    {
      SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
      if (!p) continue; // it should be GtransfoIdentity      
      p->OffsetParams(&Delta(p->Index()));
    }
}

void SimplePolyModel::FreezeErrorScales()
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& SimplePolyModel::GetTransfo(const CcdImage &Ccd) const
{
  // return GetMapping(Ccd)->Transfo(); // cannot do that
  auto p = _myMap.find(&Ccd);
  if  (p==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::GetTransfo, never heard of CcdImage "+Ccd.Name()); 
  return p->second->Transfo();  
}  

PTR(TanSipPix2RaDec) SimplePolyModel::ProduceSipWcs(const CcdImage &Ccd) const
{
  const GtransfoPoly &pix2Tp=dynamic_cast<const GtransfoPoly&>(GetTransfo(Ccd));
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(Sky2TP(Ccd));
  if (!(&pix2Tp)  || ! proj) return NULL;

  const GtransfoLin &projLinPart = proj->LinPart(); // should be the identity, but who knows? So, let us incorporate it into the pix2TP part.
  GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert())*pix2Tp;
  
  // compute a decent approximation, if higher order corrections get ignored
  GtransfoLin cdStuff = wcsPix2Tp.LinearApproximation(Ccd.ImageFrame().Center());

  // wcsPix2TP = cdStuff*sip , so
  GtransfoPoly sip = GtransfoPoly(cdStuff.invert())*wcsPix2Tp;
  Point tangentPoint( proj->TangentPoint());
  return boost::shared_ptr<TanSipPix2RaDec>(new TanSipPix2RaDec(cdStuff, tangentPoint, &sip));
}



}}}
