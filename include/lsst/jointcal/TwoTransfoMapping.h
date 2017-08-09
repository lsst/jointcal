// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H
#define LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H

#include "memory"

#include "lsst/jointcal/Mapping.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/SimplePolyMapping.h"

namespace lsst {
namespace jointcal {

//! The mapping with two transfos in a row.
class TwoTransfoMapping : public Mapping {
private:
    SimpleGtransfoMapping *_m1, *_m2;
    unsigned _nPar1, _nPar2;
    struct tmpVars  // just there to get around constness issues
    {
        Eigen::MatrixX2d h1, h2;
        Eigen::Matrix2d dt2dx;
    };

    std::unique_ptr<tmpVars> tmp;

public:
    //!
    TwoTransfoMapping(SimpleGtransfoMapping *chipMapping, SimpleGtransfoMapping *visitMapping);

    /// No copy or move: there is only ever one instance of a given model (i.e.. per ccd+visit)
    TwoTransfoMapping(TwoTransfoMapping const &) = delete;
    TwoTransfoMapping(TwoTransfoMapping &&) = delete;
    TwoTransfoMapping &operator=(TwoTransfoMapping const &) = delete;
    TwoTransfoMapping &operator=(TwoTransfoMapping &&) = delete;

    //!
    unsigned getNpar() const;

    void setMappingIndices(std::vector<unsigned> &indices) const;

    //!
    void computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint, Eigen::MatrixX2d &H) const;
    //!
    void transformPosAndErrors(FatPoint const &where, FatPoint &outPoint) const;

    /**
     * @copydoc Mapping::offsetParams
     *
     * @note  this routine is not used when fitting (the Model manages the mappings separately),
     *        but can be useful for debugging
     */
    void offsetParams(Eigen::VectorXd const &delta) {
        _m1->offsetParams(delta.segment(_m1->getIndex(), _m1->getNpar()));
        _m2->offsetParams(delta.segment(_m2->getIndex() + _m1->getNpar(), _m2->getNpar()));
    }

    //! access to transfos
    Gtransfo const &getTransfo1() const { return _m1->getTransfo(); }

    //! access to transfos
    Gtransfo const &getTransfo2() const { return _m2->getTransfo(); }

    //! Currently *not* implemented
    void positionDerivative(Point const &where, Eigen::Matrix2d &derivative, double epsilon) const;

    //! Currently not implemented
    void freezeErrorScales();

private:
    friend class ConstrainedPolyModel;
    //!
    void setWhatToFit(const bool fittingT1, const bool fittingT2);
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_TWO_TRANSFO_MAPPING_H
