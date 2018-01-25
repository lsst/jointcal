#include <memory>
#include <string>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimpleAstrometryModel.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ProjectionHandler.h"
#include "lsst/pex/exceptions.h"
#include "lsst/jointcal/Gtransfo.h"

// const int distortionDegree=3;

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.SimpleAstrometryModel");
}

namespace lsst {
namespace jointcal {

// need a way to propagate the requested degree !
SimpleAstrometryModel::SimpleAstrometryModel(CcdImageList const &ccdImageList,
                                             const std::shared_ptr<ProjectionHandler const> projectionHandler,
                                             bool initFromWcs, unsigned nNotFit, unsigned degree)
        : _sky2TP(projectionHandler)

{
    unsigned count = 0;

    for (auto i = ccdImageList.cbegin(); i != ccdImageList.cend(); ++i, ++count) {
        const CcdImage &im = **i;
        if (count < nNotFit) {
            std::unique_ptr<SimpleGtransfoMapping> id(new SimpleGtransfoMapping(GtransfoIdentity()));
            id->setIndex(-1);  // non sense, because it has no parameters
            _myMap[im.getHashKey()] = std::move(id);
        } else
        // Given how AssignIndices works, only the SimplePolyMappings
        // will actually be fitted, as nNotFit requests.
        {
            // first check that there are enough measurements for the requested polynomial degree.
            size_t nObj = im.getCatalogForFit().size();
            if (nObj == 0) {
                LOGLS_WARN(_log, "Empty catalog from image: " << im.getName());
                continue;
            }
            GtransfoPoly pol(degree);
            if (pol.getDegree() > 0)  // if not, it cannot be decreased
            {
                while (unsigned(pol.getNpar()) > 2 * nObj) {
                    LOGLS_WARN(_log, "Reducing polynomial degree from "
                                             << pol.getDegree() << ", due to too few sources (" << nObj
                                             << " vs. " << pol.getNpar() << " parameters)");
                    pol.setDegree(pol.getDegree() - 1);
                }
            }
            /* We have to center and normalize the coordinates so that
               the fit matrix is not too ill-conditionned. Basically, x
               and y in pixels are mapped to [-1,1]. When the
               transformation of SimplePolyMapping transformation is
               accessed, the combination of the normalization and the
               fitted transformation is returned, so that the trick
               remains hidden
             */
            const Frame &frame = im.getImageFrame();
            GtransfoLin shiftAndNormalize = normalizeCoordinatesTransfo(frame);
            if (initFromWcs) {
                pol = GtransfoPoly(im.getPix2TangentPlane(), frame, degree);
                pol = pol * shiftAndNormalize.invert();
            }
            _myMap[im.getHashKey()] =
                    std::unique_ptr<SimpleGtransfoMapping>(new SimplePolyMapping(shiftAndNormalize, pol));
        }
    }
}

const AstrometryMapping *SimpleAstrometryModel::getMapping(CcdImage const &ccdImage) const {
    return findMapping(ccdImage);
}

unsigned SimpleAstrometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    if (whatToFit.find("Distortions") == std::string::npos) {
        LOGLS_ERROR(_log, "AssignIndices was called and Distortions is *not* in whatToFit.");
        return 0;
    }
    unsigned index = firstIndex;
    for (auto i = _myMap.begin(); i != _myMap.end(); ++i) {
        SimplePolyMapping *p = dynamic_cast<SimplePolyMapping *>(&*(i->second));
        if (!p) continue;  // it should be GtransfoIdentity
        p->setIndex(index);
        index += p->getNpar();
    }
    return index;
}

void SimpleAstrometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

void SimpleAstrometryModel::freezeErrorTransform() {
    for (auto &i : _myMap) i.second->freezeErrorTransform();
}

const Gtransfo &SimpleAstrometryModel::getTransfo(CcdImage const &ccdImage) const {
    return dynamic_cast<const SimplePolyMapping *>(findMapping(ccdImage))->getTransfo();
}

std::shared_ptr<TanSipPix2RaDec> SimpleAstrometryModel::produceSipWcs(CcdImage const &ccdImage) const {
    const GtransfoPoly &pix2Tp = dynamic_cast<const GtransfoPoly &>(getTransfo(ccdImage));
    auto proj = std::dynamic_pointer_cast<const TanRaDec2Pix>(getSky2TP(ccdImage));
    if (!proj) return nullptr;

    const GtransfoLin &projLinPart = proj->getLinPart();  // should be the identity, but who knows? So, let us
                                                          // incorporate it into the pix2TP part.
    GtransfoPoly wcsPix2Tp = GtransfoPoly(projLinPart.invert()) * pix2Tp;

    // compute a decent approximation, if higher order corrections get ignored
    GtransfoLin cdStuff = wcsPix2Tp.linearApproximation(ccdImage.getImageFrame().getCenter());

    // wcsPix2TP = cdStuff*sip , so
    GtransfoPoly sip = GtransfoPoly(cdStuff.invert()) * wcsPix2Tp;
    Point tangentPoint(proj->getTangentPoint());
    return std::make_shared<TanSipPix2RaDec>(cdStuff, tangentPoint, &sip);
}

AstrometryMapping *SimpleAstrometryModel::findMapping(CcdImage const &ccdImage) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimpleAstrometryModel cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst