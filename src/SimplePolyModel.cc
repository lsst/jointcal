#include <string>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/SimplePolyMapping.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/Gtransfo.h"

//const int distortionDegree=3;

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.SimplePolyModel");
}

namespace lsst {
namespace jointcal {

// need a way to propagate the requested degree !
SimplePolyModel::SimplePolyModel(const CcdImageList &ccdImageList,
				 const ProjectionHandler* projectionHandler,
				 bool initFromWcs,
				 unsigned nNotFit,
				 unsigned degree) : _sky2TP(projectionHandler)

{
  // from datacards (or default)
//  unsigned degree = distortionDegree;
  unsigned count = 0;

  for (auto i=ccdImageList.cbegin(); i!= ccdImageList.cend(); ++i, ++count)
    {
      const CcdImage &im = **i;
      if (count < nNotFit)
	{
      std::unique_ptr<SimpleGtransfoMapping> id(new SimpleGtransfoMapping(GtransfoIdentity()));
	  id->setIndex(-1); // non sense, because it has no parameters
	  _myMap[&im] = std::move(id);
	}
      else
	// Given how AssignIndices works, only the SimplePolyMapping's
	// will actually be fitted, as nNotFit requests.
	{
		/* first check that there are enough measurements for the
	  requested polynomial degree */
	  size_t nObj = im.getCatalogForFit().size();
	  if (nObj == 0)
	    {
            LOGLS_WARN(_log, "Empty catalog from image: " << im.getName());
            continue;
	    }
	  GtransfoPoly pol(degree);
		if (pol.Degree() > 0) // if not, it cannot be decreased
	    while (unsigned(pol.getNpar()) > 2*nObj)
	      pol.SetDegree(pol.Degree() - 1);
	  /* We have to center and normalize the coordinates so that
	     the fit matrix is not too ill-conditionned. Basically, x
	     and y in pixels are mapped to [-1,1]. When the
	     transformation of SimplePolyMapping transformation is
	     accessed, the combination of the normalization and the
	     fitted transformation is returned, so that the trick
	     remains hidden
	   */
	  const Frame &frame  = im.getImageFrame();
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);
	  if (initFromWcs)
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


const Mapping* SimplePolyModel::getMapping(const CcdImage &ccdImageList) const
{
  mapType::const_iterator i = _myMap.find(&ccdImageList);
  if  (i==_myMap.cend()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::GetMapping, never heard of CcdImage "+ccdImageList.getName());
  return (i->second.get());
}

unsigned SimplePolyModel::assignIndices(unsigned firstIndex, const std::string &whatToFit)
{
  if (whatToFit.find("Distortions") == std::string::npos)
    {
        LOGLS_ERROR(_log, "AssignIndices was called and Distortions is *not* in whatToFit.");
        return 0;
    }
  unsigned index = firstIndex;
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    {
      SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
      if (!p) continue; // it should be GtransfoIdentity
      p->setIndex(index);
      index+= p->getNpar();
    }
  return index;
}

void SimplePolyModel::offsetParams(const Eigen::VectorXd &delta)
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    {
      SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
      if (!p) continue; // it should be GtransfoIdentity
      p->OffsetParams(&delta(p->getIndex()));
    }
}

void SimplePolyModel::freezeErrorScales()
{
  for (auto i = _myMap.begin(); i!=_myMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& SimplePolyModel::getTransfo(const CcdImage &ccdImage) const
{
  // return GetMapping(ccdImage)->Transfo(); // cannot do that
  auto p = _myMap.find(&ccdImage);
  if  (p==_myMap.end()) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,"SimplePolyModel::getTransfo, never heard of CcdImage "+ccdImage.getName());
  return p->second->Transfo();
}

std::shared_ptr<TanSipPix2RaDec> SimplePolyModel::produceSipWcs(const CcdImage &ccdImage) const
{
  const GtransfoPoly &pix2Tp=dynamic_cast<const GtransfoPoly&>(getTransfo(ccdImage));
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(sky2TP(ccdImage));
  if (!(&pix2Tp)  || ! proj) return nullptr;

  const GtransfoLin &projLinPart = proj->LinPart(); // should be the identity, but who knows? So, let us incorporate it into the pix2TP part.
  GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert())*pix2Tp;

  // compute a decent approximation, if higher order corrections get ignored
  GtransfoLin cdStuff = wcsPix2Tp.LinearApproximation(ccdImage.getImageFrame().Center());

  // wcsPix2TP = cdStuff*sip , so
  GtransfoPoly sip = GtransfoPoly(cdStuff.invert())*wcsPix2Tp;
  Point tangentPoint( proj->TangentPoint());
  return std::make_shared<TanSipPix2RaDec>(cdStuff, tangentPoint, &sip);
}

}} // end of namespaces
