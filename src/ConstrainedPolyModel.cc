#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/AstroUtils.h" // ApplyTransfo(Frame)

#include "lsst/pex/exceptions.h"
namespace pexExcept = lsst::pex::exceptions;

#include <string>
#include <iostream>

namespace lsst {
namespace jointcal {

/* This code does not contain anything involved. It just maps the
routines AstromFit needs to what is needed for this two-transfo model.
The two-transfo mappings are implemented using two one-transfo
mappings.*/

// TODO : separate the polynomial degrees for chip and visit transfos.
// TODO propagate those into python:
  static int DistortionDegree = 3;

using namespace std;

ConstrainedPolyModel::ConstrainedPolyModel(const CcdImageList &L,
					   const ProjectionHandler* ProjH,
					   bool InitFromWCS,
					   unsigned NNotFit) : _sky2TP(ProjH)

{
  // from datacards (or default)
  unsigned degree = DistortionDegree;
  unsigned count = 0;
  VisitIdType refVisit;
  // first loop to initialize all visit  and chip transfos.
  for (auto i=L.cbegin(); i!= L.cend(); ++i, ++count)
    {
      const CcdImage &im = **i;
      unsigned visit = im.getVisit();
      unsigned chip = im.getCcdId();
      auto visitp = _visitMap.find(visit);
      if (visitp == _visitMap.end())
	{
	  if (_visitMap.size() == 0)
	    {
	      // if one fits all of them, the model is degenerate.
#ifdef ROTATE_T2
# warning : hack in ConstrainedPolyModel::ConstrainedPolyModel : rotated frame
	      _visitMap[visit] = new SimpleGtransfoMapping(GtransfoLinRot(3.141927/2.), /* ToFit = */ false);
#else
	      _visitMap[visit] = std::unique_ptr<SimpleGtransfoMapping>(new SimpleGtransfoMapping(GtransfoIdentity()));
#endif
	      refVisit = visit;
	    }
	    else
#ifdef ROTATE_T2
	      {
		GtransfoPoly poly(degree);
		poly = GtransfoPoly(GtransfoLinRot(3.141927/2.))*poly;
		_visitMap[visit] = new SimplePolyMapping(GtransfoLin(), poly);
	      }
#else
	  _visitMap[visit] = std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(GtransfoLin(),
										      GtransfoPoly(degree)));
#endif
	}
      auto chipp = _chipMap.find(chip);
      if ((chipp == _chipMap.end()) && visit == refVisit )
	{
	  const Frame &frame = im.ImageFrame();

	  _tpFrame += ApplyTransfo(frame, *im.Pix2CommonTangentPlane(), LargeFrame);
	  GtransfoPoly pol(im.Pix2TangentPlane(),
			   frame,
			   degree);
	  GtransfoLin shiftAndNormalize = NormalizeCoordinatesTransfo(frame);

	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(shiftAndNormalize, pol*shiftAndNormalize.invert()));
	}
    }
  // now, second loop to set the mappings of the CCdImages
  for (auto i=L.cbegin(); i!= L.cend(); ++i, ++count)
    {
      const CcdImage &im = **i;
      unsigned visit = im.getVisit();
      unsigned chip = im.getCcdId();
      // check that the chip_indexed part was indeed assigned
      // (i.e. the reference visit was complete)
      if (_chipMap.find(chip) == _chipMap.end())
	{
	  std::cout << " WARNING: the chip " << chip << " is missing in the \
reference exposure, expect troubles" << std::endl;
	  GtransfoLin norm = NormalizeCoordinatesTransfo(im.ImageFrame());
	  _chipMap[chip] = std::unique_ptr<SimplePolyMapping>( new SimplePolyMapping(norm,
										     GtransfoPoly(degree)));
	}
      _mappings[&im] = std::unique_ptr<TwoTransfoMapping>(new TwoTransfoMapping(_chipMap[chip].get(), _visitMap[visit].get()));

    }
  cout << "INFO: ConstrainedPolyModel : we have " << _chipMap.size() << " chip mappings " << endl;
  cout << "INFO: and " << _visitMap.size() << " visit mappings " << endl;
  // DEBUG
  for (auto i=_visitMap.begin(); i != _visitMap.end(); ++i)
    cout << i-> first << endl;

}

const Mapping* ConstrainedPolyModel::GetMapping(const CcdImage &C) const
{
  mappingMapType::const_iterator i = _mappings.find(&C);
  if  (i==_mappings.end()) return nullptr;
  return (i->second.get());
}

/*! This routine decodes "DistortionsChip" and "DistortionsVisit" in
  WhatToFit. If WhatToFit contains "Distortions" and not
  Distortions<Something>, it is understood as both chips and
  visits. */
unsigned ConstrainedPolyModel::AssignIndices(unsigned FirstIndex,
					     std::string &WhatToFit)
{
  unsigned index=FirstIndex;
  if (WhatToFit.find("Distortions") == std::string::npos)
    {
      std::cout << "SimplePolyModel::AssignIndices is called and Distortions is *not*  in WhatToFit" << std::endl;
      return 0;
    }
  // if we get here "Distortions" is in WhatToFit
  _fittingChips = (WhatToFit.find("DistortionsChip") != std::string::npos);
  _fittingVisits = (WhatToFit.find("DistortionsVisit") != std::string::npos);
  // If nothing more than "Distortions" is specified, it means all:
  if ((!_fittingChips)&&(!_fittingVisits))
    {_fittingChips = _fittingVisits = true;}
  if (_fittingChips)
    for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
      {
	SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
	if (!p)
	  throw LSST_EXCEPT(pexExcept::InvalidParameterError,"ERROR: in ConstrainedPolyModel, all chip \
mappings should be SimplePolyMappings");
	p->SetIndex(index);
	index+= p->Npar();
      }
  if (_fittingVisits)
    for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
	SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
	if (!p) continue; // it should be GtransfoIdentity
	p->SetIndex(index);
	index+= p->Npar();
      }
  // Tell the mappings which derivatives they will have to fill:
  for (auto i = _mappings.begin(); i != _mappings.end() ; ++i)
    {
      i->second->SetWhatToFit(_fittingChips, _fittingVisits);
    }
  return index;
}

void ConstrainedPolyModel::OffsetParams(const Eigen::VectorXd &Delta)
{
  if (_fittingChips)
    for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->Npar()) // probably useless test
	  p->OffsetParams(&Delta(p->Index()));
      }
  if (_fittingVisits)
    for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
        auto *p = (&*(i->second));
        if (p->Npar()) // probably useless test
	  p->OffsetParams(&Delta(p->Index()));
      }
}

#if (0)
void ConstrainedPolyModel::DumpT2Transfos() const
{
   for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
      {
        auto *p = (&*(i->second));
	std::cout << "T2 for visit " << i->first << std::endl
		  << p->Transfo() << std::endl;
      }
}
#endif


void ConstrainedPolyModel::FreezeErrorScales()
{
  for (auto i = _visitMap.begin(); i!=_visitMap.end(); ++i)
    i->second->FreezeErrorScales();
  for (auto i = _chipMap.begin(); i!=_chipMap.end(); ++i)
    i->second->FreezeErrorScales();
}


const Gtransfo& ConstrainedPolyModel::GetChipTransfo(const unsigned Chip) const
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


PTR(TanSipPix2RaDec) ConstrainedPolyModel::ProduceSipWcs(const CcdImage &Ccd) const
{
  mappingMapType::const_iterator i = _mappings.find(&Ccd);
  if  (i==_mappings.end()) return nullptr;
  const TwoTransfoMapping *m = i->second.get();

  const GtransfoPoly &t1=dynamic_cast<const GtransfoPoly&>(m->T1());
  const GtransfoPoly &t2=dynamic_cast<const GtransfoPoly&>(m->T2());
  const TanRaDec2Pix *proj=dynamic_cast<const TanRaDec2Pix*>(Sky2TP(Ccd));
  if (!(&t1)  || !(&t2) || !proj) return nullptr;

  GtransfoPoly pix2Tp = t2*t1;

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
