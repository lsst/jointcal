#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/AstroUtils.h" // ApplyTransfo(Frame)

#include "lsst/pex/exceptions.h"
namespace pexExcept = lsst::pex::exceptions;

#include <string>
#include <iostream>

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.ConstrainedPolyModel");
}

namespace lsst {
namespace jointcal {

/* This code does not contain anything involved. It just maps the
routines AstrometryFit needs to what is needed for this two-transfo model.
The two-transfo mappings are implemented using two one-transfo
mappings.*/

// TODO : separate the polynomial degrees for chip and visit transfos.
// TODO propagate those into python:
  static int DistortionDegree = 3;

using namespace std;

ConstrainedPolyModel::ConstrainedPolyModel(const CcdImageList &ccdImageList,
                                           const ProjectionHandler* projectionHandler,
                                           bool initFromWCS,
                                           unsigned nNotFit) : _sky2TP(projectionHandler)

{
  // from datacards (or default)
  unsigned degree = DistortionDegree;
  // first loop to initialize all visit  and chip transfos.
  for (auto &ccdImage: ccdImageList)
    {
      const CcdImage &im = *ccdImage;
      auto visit = im.getVisit();
      auto chip = im.getCcdId();
      auto visitp = _visitMap.find(visit);
      if (visitp == _visitMap.end())
	{
	//   if (_visitMap.size() == 0)
	//     {
 //          _visitMap[visit] = std::unique_ptr<SimpleGtransfoMapping>(new SimpleGtransfoMapping(GtransfoLinScale()));
	//     }
        _visitMap[visit] = std::unique_ptr<SimpleGtransfoMapping>(new SimpleGtransfoMapping(GtransfoLin()));
	}
      auto chipp = _chipMap.find(chip);
      if (chipp == _chipMap.end())
	{
	  const Frame &frame = im.getImageFrame();

	  _tpFrame += ApplyTransfo(frame, *im.Pix2CommonTangentPlane(), LargeFrame);
	  GtransfoPoly pol(im.Pix2TangentPlane(), frame, degree);
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);

	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(shiftAndNormalize, pol*shiftAndNormalize.invert()));
	}
    }
  // now, second loop to set the mappings of the CCdImages
  for (auto &ccdImage: ccdImageList)
    {
      const CcdImage &im = *ccdImage;
      auto visit = im.getVisit();
      auto chip = im.getCcdId();
      // check that the chip_indexed part was indeed assigned
      // (i.e. the reference visit was complete)
      if (_chipMap.find(chip) == _chipMap.end())
	{
        LOGLS_WARN(_log, "Chip " << chip << " is missing in the reference exposure, expect troubles.");
	  GtransfoLin norm = NormalizeCoordinatesTransfo(im.getImageFrame());
	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>( new SimplePolyMapping(norm, GtransfoPoly(degree)));
	}
      _mappings[&im] = std::unique_ptr<TwoTransfoMapping>(new TwoTransfoMapping(_chipMap[chip].get(), _visitMap[visit].get()));

    }
  LOGLS_INFO(_log, "Constructor got " << _chipMap.size() << " chip mappings and "
             << _visitMap.size() << " visit mappings.");
  // DEBUG
  for (auto i=_visitMap.begin(); i != _visitMap.end(); ++i)
    LOGLS_DEBUG(_log, i->first);
}

const Mapping* ConstrainedPolyModel::getMapping(const CcdImage &C) const
{
  mappingMapType::const_iterator i = _mappings.find(&C);
  if  (i==_mappings.end()) return nullptr;
  return (i->second.get());
}

/*! This routine decodes "DistortionsChip" and "DistortionsVisit" in
  WhatToFit. If WhatToFit contains "Distortions" and not
  Distortions<Something>, it is understood as both chips and
  visits. */
unsigned ConstrainedPolyModel::assignIndices(unsigned FirstIndex, const std::string &WhatToFit)
{
  unsigned index=FirstIndex;
  if (WhatToFit.find("Distortions") == std::string::npos)
    {
        LOGLS_ERROR(_log, "assignIndices was called and Distortions is *not* in WhatToFit");
        return 0;
    }
  // if we get here "Distortions" is in WhatToFit
  _fittingChips = (WhatToFit.find("DistortionsChip") != std::string::npos);
  _fittingVisits = (WhatToFit.find("DistortionsVisit") != std::string::npos);
  // If nothing more than "Distortions" is specified, it means all:
  if ((!_fittingChips)&&(!_fittingVisits))
    {_fittingChips = _fittingVisits = true;}
  if (_fittingChips)
    for (auto &i: _chipMap)
      {
       i.second->setIndex(index);
       index += i.second->getNpar();
      }
  if (_fittingVisits)
    for (auto &i: _visitMap)
      {
        i.second->setIndex(index);
        index += i.second->getNpar();
      }
  // Tell the mappings which derivatives they will have to fill:
  for (auto &i: _mappings)
    {
      i.second->SetWhatToFit(_fittingChips, _fittingVisits);
    }
  return index;
}

void ConstrainedPolyModel::offsetParams(const Eigen::VectorXd &Delta)
{
  if (_fittingChips)
    for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->getNpar()) // probably useless test
	  p->OffsetParams(&Delta(p->getIndex()));
      }
  if (_fittingVisits)
    for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->getNpar()) // probably useless test
	  p->OffsetParams(&Delta(p->getIndex()));
      }
}

void ConstrainedPolyModel::freezeErrorScales()
{
  for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
    i->second->FreezeErrorScales();
  for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& ConstrainedPolyModel::getChipTransfo(const CcdIdType Chip) const
{
  auto chipp = _chipMap.find(Chip);
  if (chipp == _chipMap.end()) {
    std::stringstream errMsg;
    errMsg << "No such chipId: '" << Chip << "' found in chipMap of:  " << this;
    throw pexExcept::InvalidParameterError(errMsg.str());
  }
  return chipp->second->Transfo();
}

// Array of visits involved in the solution.
std::vector<VisitIdType> ConstrainedPolyModel::getVisits() const
{
  std::vector<VisitIdType> res;
  res.reserve(_visitMap.size());
  for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
    res.push_back(i->first);
  return res;
}

const Gtransfo& ConstrainedPolyModel::getVisitTransfo(const VisitIdType &Visit) const
{
  auto visitp = _visitMap.find(Visit);
  if (visitp == _visitMap.end()) {
    std::stringstream errMsg;
    errMsg << "No such visitId: '" << Visit << "' found in visitMap of: " << this;
    throw pexExcept::InvalidParameterError(errMsg.str());
  }
  return visitp->second->Transfo();
}


std::shared_ptr<TanSipPix2RaDec> ConstrainedPolyModel::produceSipWcs(const CcdImage &ccdImage) const
{
  const TwoTransfoMapping * mapping;
  try {
    mapping = _mappings.at(&ccdImage).get();
  }
  catch (std::out_of_range) {
    LOGLS_ERROR(_log, "CcdImage with ccd/visit " << ccdImage.getCcdId() << "/" << ccdImage.getVisit()
                << " not found in constrainedPolyModel mapping list.");
    std::ostringstream os;
    for (auto const& i: _mappings) os << i.first << ",";
    LOGLS_ERROR(_log, "Available CcdImages: " << os.str());
    return nullptr;
  }

  const GtransfoPoly &t1=dynamic_cast<const GtransfoPoly&>(mapping->T1());
  const GtransfoPoly &t2=dynamic_cast<const GtransfoPoly&>(mapping->T2());
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(sky2TP(ccdImage));
  if (!(&t1)  || !(&t2) || !proj) {
    LOGLS_ERROR(_log, "Problem with transforms of ccd/visit "
                << ccdImage.getCcdId() << "/" << ccdImage.getVisit()
                << ": T1 " << t1 << ", T2 " << t2 << ", projection " << proj);
    return nullptr;
  }

  GtransfoPoly pix2Tp = t2*t1;

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
