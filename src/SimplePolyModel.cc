#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/pex/exceptions.h"
#include <string>

#include "lsst/jointcal/Gtransfo.h"

//const int distortionDegree=3;

namespace lsst {
namespace jointcal {

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
	  id->SetIndex(-1); // non sense, because it has no parameters
	  _myMap[&im] = std::unique_ptr<SimpleGtransfoMapping>(id);
	}
      else
	// Given how AssignIndices works, only the SimplePolyMapping's
	// will actually be fitted, as NNotFit requests.
	{
		/* first check that there are enough measurements for the
	  requested polynomial degree */
	  unsigned nObj = im.getCatalogForFit().size();
	  if (nObj == 0)
	    {
	      std::cout << "WARNING: empty catalog from image : "
			<< im.getName() << std::endl;
	      continue;
	    }
	  GtransfoPoly pol(degree);
		if (pol.Degree() > 0) // if not, it cannot be decreased
	    while (unsigned(pol.Npar()) > 2*nObj)
	      pol.SetDegree(pol.Degree() - 1);
	  /* We have to center and normalize the coordinates so that
	     the fit matrix is not too ill-conditionned. Basically, x
	     and y in pixels are mapped to [-1,1]. When the
	     transformation of SimplePolyMapping transformation is
	     accessed, the combination of the normalization and the
	     fitted transformation is returned, so that the trick
	     remains hidden
	   */
	  const Frame &frame  = im.ImageFrame();
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);
	  if (InitFromWCS)
	    {
	      pol = GtransfoPoly(im.Pix2TangentPlane(),
				 frame,
				 degree);
	      pol = pol*shiftAndNormalize.invert();

	    }
	  _myMap[&im] = std::unique_ptr<SimpleGtransfoMapping>(new SimplePolyMapping(shiftAndNormalize, pol));
	}
    }
}


const Mapping* SimplePolyModel::getMapping(const CcdImage &C) const
{
  mapType::const_iterator i = _myMap.find(&C);
  if  (i==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::GetMapping, never heard of CcdImage "+C.getName());
  return (i->second.get());
}

unsigned SimplePolyModel::assignIndices(unsigned FirstIndex, std::string &WhatToFit)
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

void SimplePolyModel::offsetParams(const Eigen::VectorXd &Delta)
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    {
      SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
      if (!p) continue; // it should be GtransfoIdentity
      p->OffsetParams(&Delta(p->Index()));
    }
}

void SimplePolyModel::freezeErrorScales()
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& SimplePolyModel::GetTransfo(const CcdImage &Ccd) const
{
  // return GetMapping(Ccd)->Transfo(); // cannot do that
  auto p = _myMap.find(&Ccd);
  if  (p==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::GetTransfo, never heard of CcdImage "+Ccd.getName());
  return p->second->Transfo();
}

PTR(TanSipPix2RaDec) SimplePolyModel::ProduceSipWcs(const CcdImage &Ccd) const
{
  const GtransfoPoly &pix2Tp=dynamic_cast<const GtransfoPoly&>(GetTransfo(Ccd));
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(sky2TP(Ccd));
  if (!(&pix2Tp)  || ! proj) return nullptr;

  const GtransfoLin &projLinPart = proj->LinPart(); // should be the identity, but who knows? So, let us incorporate it into the pix2TP part.
  GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert())*pix2Tp;

  // compute a decent approximation, if higher order corrections get ignored
  GtransfoLin cdStuff = wcsPix2Tp.LinearApproximation(Ccd.ImageFrame().Center());

  // wcsPix2TP = cdStuff*sip , so
  GtransfoPoly sip = GtransfoPoly(cdStuff.invert())*wcsPix2Tp;
  Point tangentPoint( proj->TangentPoint());
  return std::shared_ptr<TanSipPix2RaDec>(new TanSipPix2RaDec(cdStuff, tangentPoint, &sip));
}

}} // end of namespaces
