#include <iostream>

#include "lsst/log/Log.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.SimplePhotometryModel");
}

namespace lsst {
namespace jointcal {

SimplePhotometryModel::SimplePhotometryModel(CcdImageList const &ccdImageList) {
    for (auto const &ccdImage : ccdImageList) {
        _myMap[ccdImage.get()] = std::unique_ptr<PhotometryMapping>(
                new PhotometryMapping(PhotometryTransfoSpatiallyInvariant()));
    }
    LOGLS_INFO(_log, "SimplePhotometryModel got " << _myMap.size() << " ccdImage mappings.");
}

unsigned SimplePhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    unsigned ipar = firstIndex;
    for (auto const &i : _myMap) {
        auto mapping = i.second.get();
        mapping->setIndex(ipar);
        ipar += mapping->getNpar();
    }
    return ipar;
}

void SimplePhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(&delta[mapping->getIndex()]);
    }
}

double SimplePhotometryModel::photomFactor(CcdImage const &ccdImage, Point const &where) const {
    auto mapping = this->findMapping(ccdImage, "photomFactor");
    return mapping->getTransfo().apply(where, 1.0);
}

void SimplePhotometryModel::getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {
    auto mapping = this->findMapping(ccdImage, "getMappingIndices");
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
}

void SimplePhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                        CcdImage const &ccdImage,
                                                        Eigen::VectorXd &derivatives) {
    // auto mapping = this->findMapping(ccdImage, "computeParameterDerivatives");
    // TODO: use mapping->computeDerivative(measuredStar)*measuredStar.getFlux() here instead.
    derivatives[0] = 1. * measuredStar.getFlux();
}

PhotometryMapping *SimplePhotometryModel::findMapping(CcdImage const &ccdImage, std::string name) const {
    auto i = _myMap.find(&ccdImage);
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimplePolyModel::" + name + ", cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
