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

#ifndef LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H
#define LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H

#include <memory>  // for unique_ptr

#include "lsst/jointcal/AstrometryMapping.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/CcdImage.h"

//! Class for a simple mapping implementing a generic AstrometryTransform
/*! It uses a template rather than a pointer so that the derived
classes can use the specifics of the transform. The class
simplePolyMapping overloads a few routines. */

namespace lsst {
namespace jointcal {

class SimpleAstrometryMapping : public AstrometryMapping {
public:
    SimpleAstrometryMapping(AstrometryTransform const &astrometryTransform, bool toBeFit = true)
            : toBeFit(toBeFit),
              transform(astrometryTransform.clone()),
              errorProp(transform),
              lin(new AstrometryTransformLinear) {
        // in this order:
        // take a copy of the input transform,
        // assign the transformation used to propagate errors to the transform itself
        // reserve some memory space to compute the derivatives (efficiency).
    }

    /// No copy or move: there is only ever one instance of a given mapping (i.e.. per ccd+visit)
    SimpleAstrometryMapping(SimpleAstrometryMapping const &) = delete;
    SimpleAstrometryMapping(SimpleAstrometryMapping &&) = delete;
    SimpleAstrometryMapping &operator=(SimpleAstrometryMapping const &) = delete;
    SimpleAstrometryMapping &operator=(SimpleAstrometryMapping &&) = delete;

    virtual void freezeErrorTransform() {
        // from there on, updating the transform does not change the errors.
        errorProp = transform->clone();
    }

    // interface Mapping functions:

    //!
    unsigned getNpar() const {
        if (toBeFit)
            return transform->getNpar();
        else
            return 0;
    }

    //!
    void getMappingIndices(std::vector<unsigned> &indices) const {
        if (indices.size() < getNpar()) indices.resize(getNpar());
        for (unsigned k = 0; k < getNpar(); ++k) indices[k] = index + k;
    }

    //!
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const {
        transform->transformPosAndErrors(where, outPoint);
        FatPoint tmp;
        errorProp->transformPosAndErrors(where, tmp);
        outPoint.vx = tmp.vx;
        outPoint.vy = tmp.vy;
        outPoint.vxy = tmp.vxy;
    }

    //!
    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const {
        errorProp->computeDerivative(where, *lin, epsilon);
        derivative(0, 0) = lin->coeff(1, 0, 0);
        //
        /* This does not work : it was proved by rotating the frame
           see the compilation switch ROTATE_T2 in constrainedAstrometryModel.cc
        derivative(1,0) = lin->coeff(1,0,1);
        derivative(0,1) = lin->coeff(0,1,0);
        */
        derivative(1, 0) = lin->coeff(0, 1, 0);
        derivative(0, 1) = lin->coeff(1, 0, 1);
        derivative(1, 1) = lin->coeff(0, 1, 1);
    }

    //!
    void offsetParams(Eigen::VectorXd const &delta) {
        if (toBeFit) transform->offsetParams(delta);
    }

    //! position of the parameters within the grand fitting scheme
    unsigned getIndex() const { return index; }

    //!
    void setIndex(unsigned i) { index = i; }

    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const {
        transformPosAndErrors(where, outPoint);
        transform->paramDerivatives(where, &H(0, 0), &H(0, 1));
    }

    //! Access to the (fitted) transform
    virtual AstrometryTransform const &getTransform() const { return *transform; }

    /// Get whether this mapping is fit as part of a Model.
    bool getToBeFit() const { return toBeFit; }
    /// Set whether this Mapping is to be fit as part of a Model.
    void setToBeFit(bool value) { toBeFit = value; }

protected:
    // Whether this Mapping is fit as part of a Model.
    bool toBeFit;
    unsigned index;
    /* inheritance may also work. Perhaps with some trouble because
       some routines in Mapping and AstrometryTransform have the same name */
    std::shared_ptr<AstrometryTransform> transform;

    std::shared_ptr<AstrometryTransform> errorProp;
    /* to avoid allocation at every call of positionDerivative.
       use a pointer for constness */
    std::unique_ptr<AstrometryTransformLinear> lin;
};

//! Mapping implementation for a polynomial transformation.
class SimplePolyMapping : public SimpleAstrometryMapping {
public:
    ~SimplePolyMapping() {}

    // ! contructor.
    /*! The transformation will be initialized to transform, so that the effective transformation
      reads transform*CenterAndScale */
    SimplePolyMapping(AstrometryTransformLinear const &CenterAndScale,
                      AstrometryTransformPolynomial const &transform)
            : SimpleAstrometryMapping(transform), _centerAndScale(CenterAndScale) {
        // We assume that the initialization was done properly, for example that
        // transform = pixToTangentPlane*CenterAndScale.inverted(), so we do not touch transform.
        /* store the (spatial) derivative of _centerAndScale. For the extra
           diagonal terms, just copied the ones in positionDerivatives */
        preDer(0, 0) = _centerAndScale.coeff(1, 0, 0);
        preDer(1, 0) = _centerAndScale.coeff(0, 1, 0);
        preDer(0, 1) = _centerAndScale.coeff(1, 0, 1);
        preDer(1, 1) = _centerAndScale.coeff(0, 1, 1);

        // check of matrix indexing (once for all)
        MatrixX2d H(3, 2);
        assert((&H(1, 0) - &H(0, 0)) == 1);
    }

    /// No copy or move: there is only ever one instance of a given mapping (i.e.. per ccd+visit)
    SimplePolyMapping(SimplePolyMapping const &) = delete;
    SimplePolyMapping(SimplePolyMapping &&) = delete;
    SimplePolyMapping &operator=(SimplePolyMapping const &) = delete;
    SimplePolyMapping &operator=(SimplePolyMapping &&) = delete;

    /* The SimpleAstrometryMapping version does not account for the
       _centerAndScale transform */

    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const {
        Point tmp = _centerAndScale.apply(where);
        errorProp->computeDerivative(tmp, *lin, epsilon);
        derivative(0, 0) = lin->coeff(1, 0, 0);
        //
        /* This does not work : it was proved by rotating the frame
           see the compilation switch ROTATE_T2 in constrainedAstrometryModel.cc
        derivative(1,0) = lin->coeff(1,0,1);
        derivative(0,1) = lin->coeff(0,1,0);
        */
        derivative(1, 0) = lin->coeff(0, 1, 0);
        derivative(0, 1) = lin->coeff(1, 0, 1);
        derivative(1, 1) = lin->coeff(0, 1, 1);
        derivative = preDer * derivative;
    }

    //! Calls the transforms and implements the centering and scaling of coordinates
    /* We should put the computation of error propagation and
       parameter derivatives into the same AstrometryTransform routine because
       it could be significantly faster */
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
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

    //! Implements as well the centering and scaling of coordinates
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const {
        FatPoint mid;
        _centerAndScale.transformPosAndErrors(where, mid);
        transform->transformPosAndErrors(mid, outPoint);
        FatPoint tmp;
        errorProp->transformPosAndErrors(mid, tmp);
        outPoint.vx = tmp.vx;
        outPoint.vy = tmp.vy;
        outPoint.vxy = tmp.vxy;
    }

    //! Access to the (fitted) transform
    AstrometryTransform const &getTransform() const {
        // Cannot fail given the contructor:
        const AstrometryTransformPolynomial *fittedPoly =
                dynamic_cast<const AstrometryTransformPolynomial *>(&(*transform));
        actualResult = (*fittedPoly) * _centerAndScale;
        return actualResult;
    }

private:
    /* to better condition the 2nd derivative matrix, the
    transformed coordinates are mapped (roughly) on [-1,1].
    We need both the transform and its derivative. */
    AstrometryTransformLinear _centerAndScale;
    Eigen::Matrix2d preDer;

    /* Where we store the combination. */
    mutable AstrometryTransformPolynomial actualResult;
};

#ifdef STORAGE
/*! "do nothing" mapping. The Ccdimage's that "use" this one impose the
  coordinate system */
class SimpleIdentityMapping : public SimpleAstrometryMapping<AstrometryTransformIdentity> {
public:
    //! nothing to do.
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const {
        outPoint = where;
    }
};
#endif
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_ASTROMETRY_MAPPING_H
