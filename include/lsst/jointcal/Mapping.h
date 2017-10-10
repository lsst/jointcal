// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_MAPPING_H
#define LSST_JOINTCAL_MAPPING_H

#include <vector>
#include "lsst/jointcal/Eigenstuff.h"

namespace lsst {
namespace jointcal {

class FatPoint;
class Point;

//! virtual class needed in the abstraction of the distortion model
class Mapping {
public:
    //! Number of parameters in total
    virtual unsigned getNpar() const = 0;

    /// Sets how this set of parameters (of length Npar()) map into the "grand" fit
    /// Expects that indices has enough space reserved.
    virtual void getMappingIndices(std::vector<unsigned> &indices) const = 0;

    //! Actually applies the mapping and evaluates the derivatives w.r.t the fitted parameters.
    /*! This is grouped into a single call because for most models,
        evaluating the derivatives w.r.T parameters is not much longer
        than just transforming */
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const = 0;
    //! The same as above but without the parameter derivatives (used to evaluate chi^2)
    virtual void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const = 0;

    //! Remember the error scale and freeze it
    //  virtual void freezeErrorScales() = 0;

    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    //! The derivative w.r.t. position
    virtual void positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                    double epsilon) const = 0;

    //!
    virtual ~Mapping(){};
};
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_MAPPING_H
