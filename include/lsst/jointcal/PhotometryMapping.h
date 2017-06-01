// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
#define LSST_JOINTCAL_PHOTOMETRY_MAPPING_H

#include <memory>

#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

class Point;

class PhotometryMapping {
public:
    PhotometryMapping(PhotometryTransfo const &_transfo) : index(-1), transfo(_transfo.clone()) {}

    /// Number of total parameters in this mapping
    unsigned getNpar() const { return transfo->getNpar(); }

    /*
     * Sets how this set of parameters (of length Npar()) map into the "grand" fit.
     * Expects that indices has enough space reserved.
     */
    void setMappingIndices(std::vector<unsigned> &indices) const {
        indices.reserve(getNpar());
        for (unsigned k = 0; k < getNpar(); ++k) indices[k] = index + k;
    }

    /**
     * Applies the mapping and evaluates the derivatives with respect to the fitted parameters.
     *
     * This is grouped into a single call because for most models,
     * evaluating the derivatives w.r.t parameters is not much longer than just transforming.
     */
    void computeTransformAndDerivatives(const Point &where, double &out, Eigen::MatrixX2d &H) const;

    //! The same as above but without the parameter derivatives (used to evaluate chi^2)
    void transformPosAndErrors(const Point &where, double &out) const;

    void offsetParams(const double *delta) { transfo->offsetParams(delta); }

    /// Get the index of this mapping in the grand fit.
    unsigned getIndex() { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(unsigned i) { index = i; }

    PhotometryTransfo const &getTransfo() { return *(transfo.get()); }

protected:
    // Start index of this mapping in the "grand" fit
    unsigned index;

    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransfo> transfo;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
