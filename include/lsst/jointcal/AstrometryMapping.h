// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef LSST_JOINTCAL_ASTROMETRY_MAPPING_H
#define LSST_JOINTCAL_ASTROMETRY_MAPPING_H

#include <iostream>
#include <vector>
#include "lsst/jointcal/Eigenstuff.h"

namespace lsst {
namespace jointcal {

class FatPoint;
class Point;

//! virtual class needed in the abstraction of the distortion model
class AstrometryMapping {
public:
    //! Number of parameters in total
    virtual std::size_t getNpar() const = 0;

    /// Sets how this set of parameters (of length Npar()) map into the "grand" fit
    /// Expects that indices has enough space reserved.
    virtual void getMappingIndices(IndexVector &indices) const = 0;

    //! Actually applies the AstrometryMapping and evaluates the derivatives w.r.t the fitted parameters.
    /*! This is grouped into a single call because for most models,
        evaluating the derivatives w.r.T parameters is not much longer
        than just transforming */
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const = 0;
    //! The same as above but without the parameter derivatives (used to evaluate chi^2)
    virtual void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const = 0;

    //! Remember the error scale and freeze it
    //  virtual void freezeErrorTransform() = 0;

    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    //! The derivative w.r.t. position
    virtual void positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                    double epsilon) const = 0;

    /**
     * Print a string representation of the contents of this mapping, for debugging.
     *
     * This string representation can be very verbose, as it contains all of the parameters
     * of all of the transforms in this mapping.
     */
    virtual void print(std::ostream &out) const = 0;

    //!
    virtual ~AstrometryMapping(){};
};

inline std::ostream &operator<<(std::ostream &stream, AstrometryMapping const &mapping) {
    mapping.print(stream);
    return stream;
}

}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_ASTROMETRY_MAPPING_H
