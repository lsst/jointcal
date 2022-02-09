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

#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/SimpleAstrometryMapping.h"

namespace lsst {
namespace jointcal {

void SimpleAstrometryMapping::getMappingIndices(IndexVector &indices) const {
    if (indices.size() < getNpar()) {
        indices.resize(getNpar());
    }
    for (std::size_t k = 0; k < getNpar(); ++k) {
        indices[k] = index + k;
    }
}

void SimpleAstrometryMapping::transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const {
    transform->transformPosAndErrors(where, outPoint);
    FatPoint tmp;
    errorProp->transformPosAndErrors(where, tmp);
    outPoint.vx = tmp.vx;
    outPoint.vy = tmp.vy;
    outPoint.vxy = tmp.vxy;
}

void SimpleAstrometryMapping::positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                                 double epsilon) const {
    errorProp->computeDerivative(where, *lin, epsilon);
    derivative(0, 0) = lin->getCoefficient(1, 0, 0);
    //
    /* This does not work : it was proved by rotating the frame
       see the compilation switch ROTATE_T2 in constrainedAstrometryModel.cc
    derivative(1,0) = lin->getCoefficient(1,0,1);
    derivative(0,1) = lin->getCoefficient(0,1,0);
    */
    derivative(1, 0) = lin->getCoefficient(0, 1, 0);
    derivative(0, 1) = lin->getCoefficient(1, 0, 1);
    derivative(1, 1) = lin->getCoefficient(0, 1, 1);
}

void SimpleAstrometryMapping::computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                             Eigen::MatrixX2d &H) const {
    transformPosAndErrors(where, outPoint);
    transform->paramDerivatives(where, &H(0, 0), &H(0, 1));
}

void SimpleAstrometryMapping::print(std::ostream &out) const { out << *transform; }

SimplePolyMapping::SimplePolyMapping(AstrometryTransformLinear const &CenterAndScale,
                                     AstrometryTransformPolynomial const &transform)
        : SimpleAstrometryMapping(transform), _centerAndScale(CenterAndScale) {
    // We assume that the initialization was done properly, for example that
    // transform = pixToTangentPlane*CenterAndScale.inverted(), so we do not touch transform.
    /* store the (spatial) derivative of _centerAndScale. For the extra
       diagonal terms, just copied the ones in positionDerivatives */
    preDer(0, 0) = _centerAndScale.getCoefficient(1, 0, 0);
    preDer(1, 0) = _centerAndScale.getCoefficient(0, 1, 0);
    preDer(0, 1) = _centerAndScale.getCoefficient(1, 0, 1);
    preDer(1, 1) = _centerAndScale.getCoefficient(0, 1, 1);

    // check of matrix indexing (once for all)
    MatrixX2d H(3, 2);
    assert((&H(1, 0) - &H(0, 0)) == 1);
}

void SimplePolyMapping::positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                           double epsilon) const {
    Point tmp = _centerAndScale.apply(where);
    errorProp->computeDerivative(tmp, *lin, epsilon);
    derivative(0, 0) = lin->getCoefficient(1, 0, 0);
    //
    /* This does not work : it was proved by rotating the frame
       see the compilation switch ROTATE_T2 in constrainedAstrometryModel.cc
    derivative(1,0) = lin->getCoefficient(1,0,1);
    derivative(0,1) = lin->getCoefficient(0,1,0);
    */
    derivative(1, 0) = lin->getCoefficient(0, 1, 0);
    derivative(0, 1) = lin->getCoefficient(1, 0, 1);
    derivative(1, 1) = lin->getCoefficient(0, 1, 1);
    derivative = preDer * derivative;
}

void SimplePolyMapping::computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                       Eigen::MatrixX2d &H) const {
    FatPoint mid;
    _centerAndScale.transformPosAndErrors(where, mid);
    transform->transformPosAndErrors(mid, outPoint);
    FatPoint tmp;
    errorProp->transformPosAndErrors(mid, tmp);
    outPoint.vx = tmp.vx;
    outPoint.vy = tmp.vy;
    outPoint.vxy = tmp.vxy;
    transform->paramDerivatives(mid, &H(0, 0), &H(0, 1));
}

void SimplePolyMapping::transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const {
    FatPoint mid;
    _centerAndScale.transformPosAndErrors(where, mid);
    transform->transformPosAndErrors(mid, outPoint);
    FatPoint tmp;
    errorProp->transformPosAndErrors(mid, tmp);
    outPoint.vx = tmp.vx;
    outPoint.vy = tmp.vy;
    outPoint.vxy = tmp.vxy;
}

AstrometryTransform const &SimplePolyMapping::getTransform() const {
    // Cannot fail given the contructor:
    const auto *fittedPoly =
            dynamic_cast<const AstrometryTransformPolynomial *>(&(*transform));
    actualResult = (*fittedPoly) * _centerAndScale;
    return actualResult;
}

}  // namespace jointcal
}  // namespace lsst
