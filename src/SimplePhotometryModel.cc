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
        _myMap[ccdImage.get()] =
                std::unique_ptr<PhotometryMapping>(new PhotometryMapping(ConstantPhotometryTransfo()));
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
        mapping->offsetParams(&delta(mapping->getIndex()));
    }
}

const PhotometryTransfo &SimplePhotometryModel::getTransfo(CcdImage const &ccdImage) const {
    auto mapping = this->findMapping(ccdImage, "getTransfo");
    return mapping->getTransfo();
}

double SimplePhotometryModel::photomFactor(CcdImage const &ccdImage, Point const &where) const {
    auto mapping = this->findMapping(ccdImage, "photomFactor");
    return mapping->getTransfo().apply(where, 1.0);
}

void SimplePhotometryModel::getIndicesAndDerivatives(MeasuredStar const &measuredStar,
                                                     CcdImage const &ccdImage, std::vector<unsigned> &indices,
                                                     Eigen::VectorXd &derivative) {
    auto mapping = this->findMapping(ccdImage, "getIndicesAndDerivatives");
    indices.resize(1);
    indices[0] = mapping->getIndex();
    derivative[0] = 1;
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
