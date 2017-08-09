#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstrometryModel.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/jointcal/AstroUtils.h"  // applyTransfo(Frame)

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

ConstrainedPolyModel::ConstrainedPolyModel(CcdImageList const &ccdImageList,
                                           ProjectionHandler const *projectionHandler, bool initFromWCS,
                                           unsigned nNotFit)
        : _sky2TP(projectionHandler)

{
    // from datacards (or default)
    unsigned degree = DistortionDegree;
    // first loop to initialize all visit  and chip transfos.
    for (auto &ccdImage : ccdImageList) {
        const CcdImage &im = *ccdImage;
        auto visit = im.getVisit();
        auto chip = im.getCcdId();
        auto visitp = _visitMap.find(visit);
        if (visitp == _visitMap.end()) {
            if (_visitMap.size() == 0) {
                _visitMap[visit] =
                        std::unique_ptr<SimpleGtransfoMapping>(new SimpleGtransfoMapping(GtransfoIdentity()));
            } else {
                _visitMap[visit] = std::unique_ptr<SimpleGtransfoMapping>(
                        new SimplePolyMapping(GtransfoLin(), GtransfoPoly(degree)));
            }
        }
        auto chipp = _chipMap.find(chip);
        if (chipp == _chipMap.end()) {
            const Frame &frame = im.getImageFrame();

            _tpFrame += applyTransfo(frame, *im.getPix2CommonTangentPlane(), LargeFrame);
            GtransfoPoly pol(im.getPix2TangentPlane(), frame, degree);
            GtransfoLin shiftAndNormalize = normalizeCoordinatesTransfo(frame);

            _chipMap[chip] = std::unique_ptr<SimplePolyMapping>(
                    new SimplePolyMapping(shiftAndNormalize, pol * shiftAndNormalize.invert()));
        }
    }
    // now, second loop to set the mappings of the CCdImages
    for (auto &ccdImage : ccdImageList) {
        const CcdImage &im = *ccdImage;
        auto visit = im.getVisit();
        auto chip = im.getCcdId();
        // check that the chip_indexed part was indeed assigned
        // (i.e. the reference visit was complete)
        if (_chipMap.find(chip) == _chipMap.end()) {
            LOGLS_WARN(_log, "Chip " << chip << " is missing in the reference exposure, expect troubles.");
            GtransfoLin norm = normalizeCoordinatesTransfo(im.getImageFrame());
            _chipMap[chip] =
                    std::unique_ptr<SimplePolyMapping>(new SimplePolyMapping(norm, GtransfoPoly(degree)));
        }
        _mappings[&im] = std::unique_ptr<TwoTransfoMapping>(
                new TwoTransfoMapping(_chipMap[chip].get(), _visitMap[visit].get()));
    }
    LOGLS_INFO(_log, "Constructor got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                                        << " visit mappings.");
    // DEBUG
    for (auto i = _visitMap.begin(); i != _visitMap.end(); ++i) LOGLS_DEBUG(_log, i->first);
}

const Mapping *ConstrainedPolyModel::getMapping(CcdImage const &ccdImage) const {
    mappingMapType::const_iterator i = _mappings.find(&ccdImage);
    if (i == _mappings.end()) return nullptr;
    return (i->second.get());
}

/*! This routine decodes "DistortionsChip" and "DistortionsVisit" in
  whatToFit. If whatToFit contains "Distortions" and not
  Distortions<Something>, it is understood as both chips and
  visits. */
unsigned ConstrainedPolyModel::assignIndices(unsigned firstIndex, std::string const &whatToFit) {
    unsigned index = firstIndex;
    if (whatToFit.find("Distortions") == std::string::npos) {
        LOGLS_ERROR(_log, "assignIndices was called and Distortions is *not* in whatToFit");
        return 0;
    }
    // if we get here "Distortions" is in whatToFit
    _fittingChips = (whatToFit.find("DistortionsChip") != std::string::npos);
    _fittingVisits = (whatToFit.find("DistortionsVisit") != std::string::npos);
    // If nothing more than "Distortions" is specified, it means all:
    if ((!_fittingChips) && (!_fittingVisits)) {
        _fittingChips = _fittingVisits = true;
    }
    if (_fittingChips)
        for (auto &i : _chipMap) {
            i.second->setIndex(index);
            index += i.second->getNpar();
        }
    if (_fittingVisits)
        for (auto &i : _visitMap) {
            i.second->setIndex(index);
            index += i.second->getNpar();
        }
    // Tell the mappings which derivatives they will have to fill:
    for (auto &i : _mappings) {
        i.second->setWhatToFit(_fittingChips, _fittingVisits);
    }
    return index;
}

void ConstrainedPolyModel::offsetParams(Eigen::VectorXd const &delta) {
    if (_fittingChips)
        for (auto &i : _chipMap) {
            auto mapping = i.second.get();
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
    if (_fittingVisits)
        for (auto &i : _visitMap) {
            auto mapping = i.second.get();
            mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
        }
}

void ConstrainedPolyModel::freezeErrorScales() {
    for (auto i = _visitMap.begin(); i != _visitMap.end(); ++i) i->second->freezeErrorScales();
    for (auto i = _chipMap.begin(); i != _chipMap.end(); ++i) i->second->freezeErrorScales();
}

const Gtransfo &ConstrainedPolyModel::getChipTransfo(CcdIdType const chip) const {
    auto chipp = _chipMap.find(chip);
    if (chipp == _chipMap.end()) {
        std::stringstream errMsg;
        errMsg << "No such chipId: '" << chip << "' found in chipMap of:  " << this;
        throw pexExcept::InvalidParameterError(errMsg.str());
    }
    return chipp->second->getTransfo();
}

// Array of visits involved in the solution.
std::vector<VisitIdType> ConstrainedPolyModel::getVisits() const {
    std::vector<VisitIdType> res;
    res.reserve(_visitMap.size());
    for (auto i = _visitMap.begin(); i != _visitMap.end(); ++i) res.push_back(i->first);
    return res;
}

const Gtransfo &ConstrainedPolyModel::getVisitTransfo(VisitIdType const &visit) const {
    auto visitp = _visitMap.find(visit);
    if (visitp == _visitMap.end()) {
        std::stringstream errMsg;
        errMsg << "No such visitId: '" << visit << "' found in visitMap of: " << this;
        throw pexExcept::InvalidParameterError(errMsg.str());
    }
    return visitp->second->getTransfo();
}

std::shared_ptr<TanSipPix2RaDec> ConstrainedPolyModel::produceSipWcs(CcdImage const &ccdImage) const {
    const TwoTransfoMapping *mapping;
    try {
        mapping = _mappings.at(&ccdImage).get();
    } catch (std::out_of_range &) {
        LOGLS_ERROR(_log, "CcdImage with ccd/visit " << ccdImage.getCcdId() << "/" << ccdImage.getVisit()
                                                     << " not found in constrainedPolyModel mapping list.");
        std::ostringstream os;
        for (auto const &i : _mappings) os << i.first << ",";
        LOGLS_ERROR(_log, "Available CcdImages: " << os.str());
        return nullptr;
    }

    GtransfoPoly pix2Tp;
    const GtransfoPoly &t1 = dynamic_cast<const GtransfoPoly &>(mapping->getTransfo1());
    // TODO: This line produces a warning on clang (t1 is always valid: a failed dynamic_cast of a reference
    // raises bad_cast instead of returning nullptr like a failed pointer cast), but I'll deal with it as
    // part of DM-10524 (hopefully removing the necessity of the casts).
    if (!(&t1)) {
        LOGLS_ERROR(_log, "Problem with transform 1 of ccd/visit " << ccdImage.getCcdId() << "/"
                                                                   << ccdImage.getVisit() << ": T1 "
                                                                   << mapping->getTransfo1());
        return nullptr;
    }
    // NOTE: we currently expect T2 to be an identity for the first visit, so we have to treat it separately.
    // TODO: We are aware that this is a hack, but it will be fixed as part of DM-10524.
    try {
        const GtransfoIdentity &t2 = dynamic_cast<const GtransfoIdentity &>(mapping->getTransfo2());
        pix2Tp = t1;
    } catch (std::bad_cast &) {
        try {
            const GtransfoPoly &t2_poly = dynamic_cast<const GtransfoPoly &>(mapping->getTransfo2());
            pix2Tp = t1 * t2_poly;
        } catch (std::bad_cast &) {
            LOGLS_ERROR(_log, "Problem with transform 2 of ccd/visit " << ccdImage.getCcdId() << "/"
                                                                       << ccdImage.getVisit() << ": T2 "
                                                                       << mapping->getTransfo2());
            return nullptr;
        }
    }
    const TanRaDec2Pix *proj = dynamic_cast<const TanRaDec2Pix *>(getSky2TP(ccdImage));
    if (!proj) {
        LOGLS_ERROR(_log, "Problem with projection of ccd/visit " << ccdImage.getCcdId() << "/"
                                                                  << ccdImage.getVisit() << ": projection "
                                                                  << getSky2TP(ccdImage));
        return nullptr;
    }

    // should be the identity, but who knows? So, let us incorporate it into the pix2TP part.
    const GtransfoLin &projLinPart = proj->getLinPart();
    GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert()) * pix2Tp;

    // compute a decent approximation, if higher order corrections get ignored
    GtransfoLin cdStuff = wcsPix2Tp.linearApproximation(ccdImage.getImageFrame().getCenter());

    // wcsPix2TP = cdStuff*sip , so
    GtransfoPoly sip = GtransfoPoly(cdStuff.invert()) * wcsPix2Tp;
    Point tangentPoint(proj->getTangentPoint());
    return std::make_shared<TanSipPix2RaDec>(cdStuff, tangentPoint, &sip);
}
}  // namespace jointcal
}  // namespace lsst
