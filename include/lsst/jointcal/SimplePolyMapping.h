// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_POLY_MAPPING_H
#define LSST_JOINTCAL_SIMPLE_POLY_MAPPING_H

#include <memory>  // for unique_ptr

#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/CcdImage.h"

//! Class for a simple mapping implementing a generic Gtransfo
/*! It uses a template rather than a pointer so that the derived
classes can use the specifics of the transfo. The class
simplePolyMapping overloads a few routines. */

namespace lsst {
namespace jointcal {

class SimpleGtransfoMapping : public Mapping {
protected:
    bool toFit;
    unsigned index;
    /* inheritance may also work. Perhaps with some trouble because
       some routines in Mapping and Gtransfo have the same name */
    std::shared_ptr<Gtransfo> transfo;

    std::shared_ptr<Gtransfo> errorProp;
    /* to avoid allocation at every call of positionDerivative.
       use a pointer for constness */
    std::unique_ptr<GtransfoLin> lin;

public:
    SimpleGtransfoMapping(Gtransfo const &gtransfo, bool toFit = true)
            : toFit(toFit), transfo(gtransfo.clone()), errorProp(transfo), lin(new GtransfoLin) {
        // in this order:
        // take a copy of the input transfo,
        // assign the transformation used to propagate errors to the transfo itself
        // reserve some memory space to compute the derivatives (efficiency).
    }

    /// No copy or move: there is only ever one instance of a given mapping (i.e.. per ccd+visit)
    SimpleGtransfoMapping(SimpleGtransfoMapping const &) = delete;
    SimpleGtransfoMapping(SimpleGtransfoMapping &&) = delete;
    SimpleGtransfoMapping &operator=(SimpleGtransfoMapping const &) = delete;
    SimpleGtransfoMapping &operator=(SimpleGtransfoMapping &&) = delete;

    virtual void freezeErrorScales() {
        // from there on, updating the transfo does not change the errors.
        errorProp = transfo->clone();
    }

    // interface Mapping functions:

    //!
    unsigned getNpar() const {
        if (toFit)
            return transfo->getNpar();
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
        transfo->transformPosAndErrors(where, outPoint);
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
           see the compilation switch ROTATE_T2 in constrainedpolymodel.cc
        derivative(1,0) = lin->coeff(1,0,1);
        derivative(0,1) = lin->coeff(0,1,0);
        */
        derivative(1, 0) = lin->coeff(0, 1, 0);
        derivative(0, 1) = lin->coeff(1, 0, 1);
        derivative(1, 1) = lin->coeff(0, 1, 1);
    }

    //!
    void offsetParams(Eigen::VectorXd const &delta) { transfo->offsetParams(delta); }

    //! position of the parameters within the grand fitting scheme
    unsigned getIndex() const { return index; }

    //!
    void setIndex(unsigned i) { index = i; }

    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const {
        transformPosAndErrors(where, outPoint);
        transfo->paramDerivatives(where, &H(0, 0), &H(0, 1));
    }

    //! Access to the (fitted) transfo
    virtual Gtransfo const &getTransfo() const { return *transfo; }
};

//! Mapping implementation for a polynomial transformation.
class SimplePolyMapping : public SimpleGtransfoMapping {
    /* to better condition the 2nd derivative matrix, the
    transformed coordinates are mapped (roughly) on [-1,1].
    We need both the transform and its derivative. */
    GtransfoLin _centerAndScale;
    Eigen::Matrix2d preDer;

    /* Where we store the combination. */
    mutable GtransfoPoly actualResult;

public:
    ~SimplePolyMapping() {}

    // ! contructor.
    /*! The transformation will be initialized to gtransfo, so that the effective transformation
      reads gtransfo*CenterAndScale */
    SimplePolyMapping(GtransfoLin const &CenterAndScale, GtransfoPoly const &gtransfo)
            : SimpleGtransfoMapping(gtransfo), _centerAndScale(CenterAndScale) {
        // We assume that the initialization was done properly, for example that
        // gtransfo = pix2TP*CenterAndScale.invert(), so we do not touch transfo.
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

    /* The SimpleGtransfoMapping version does not account for the
       _centerAndScale transfo */

    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const {
        Point tmp = _centerAndScale.apply(where);
        errorProp->computeDerivative(tmp, *lin, epsilon);
        derivative(0, 0) = lin->coeff(1, 0, 0);
        //
        /* This does not work : it was proved by rotating the frame
           see the compilation switch ROTATE_T2 in constrainedpolymodel.cc
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
       parameter derivatives into the same Gtransfo routine because
       it could be significantly faster */
    virtual void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                Eigen::MatrixX2d &H) const {
        FatPoint mid;
        _centerAndScale.transformPosAndErrors(where, mid);
        transfo->transformPosAndErrors(mid, outPoint);
        FatPoint tmp;
        errorProp->transformPosAndErrors(mid, tmp);
        outPoint.vx = tmp.vx;
        outPoint.vy = tmp.vy;
        outPoint.vxy = tmp.vxy;
        transfo->paramDerivatives(mid, &H(0, 0), &H(0, 1));
    }

    //! Implements as well the centering and scaling of coordinates
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const {
        FatPoint mid;
        _centerAndScale.transformPosAndErrors(where, mid);
        transfo->transformPosAndErrors(mid, outPoint);
        FatPoint tmp;
        errorProp->transformPosAndErrors(mid, tmp);
        outPoint.vx = tmp.vx;
        outPoint.vy = tmp.vy;
        outPoint.vxy = tmp.vxy;
    }

    //! Access to the (fitted) transfo
    Gtransfo const &getTransfo() const {
        // Cannot fail given the contructor:
        const GtransfoPoly *fittedPoly = dynamic_cast<const GtransfoPoly *>(&(*transfo));
        actualResult = (*fittedPoly) * _centerAndScale;
        return actualResult;
    }
};

#ifdef STORAGE
/*! "do nothing" mapping. The Ccdimage's that "use" this one impose the
  coordinate system */
class SimpleIdentityMapping : public SimpleGtransfoMapping<GtransfoIdentity> {
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

#endif  // LSST_JOINTCAL_SIMPLE_POLY_MAPPING_H
