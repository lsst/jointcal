#include "lsst/jointcal/TwoTransfoMapping.h"
#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

namespace lsst {
namespace jointcal {

TwoTransfoMapping::TwoTransfoMapping(SimpleGtransfoMapping *mapping1, SimpleGtransfoMapping *mapping2)
        : _m1(mapping1), _m2(mapping2) {
    /* Allocate the record of temporary variables, so that they are not
       allocated at every call. This is hidden behind a pointer in order
       to be allowed to alter them in a const routine. */
    tmp = std::unique_ptr<tmpVars>(new tmpVars);
    setWhatToFit(true, true);
}

unsigned TwoTransfoMapping::getNpar() const { return _nPar1 + _nPar2; }

void TwoTransfoMapping::getMappingIndices(std::vector<unsigned> &indices) const {
    unsigned npar = getNpar();
    if (indices.size() < npar) indices.resize(npar);
    // in case we are only fitting one of the two transfos
    if (_nPar1)
        _m1->getMappingIndices(indices);
    else if (_nPar2) {
        _m2->getMappingIndices(indices);
        return;
    }
    // if we get here we are fitting both
    // there is probably a more elegant way to feed a subpart of a std::vector
    std::vector<unsigned> ind2(_nPar2);
    _m2->getMappingIndices(ind2);
    for (unsigned k = 0; k < _nPar2; ++k) indices.at(k + _nPar1) = ind2.at(k);
}

void TwoTransfoMapping::computeTransformAndDerivatives(FatPoint const &where, FatPoint &outPoint,
                                                       Eigen::MatrixX2d &H) const {
    // not true in general. Will crash if H is too small.
    //  assert(H.cols()==Npar());

    FatPoint pMid;
    // don't need errors there but no Mapping::Transform() routine.

    if (_nPar1) {
        _m1->computeTransformAndDerivatives(where, pMid, tmp->h1);
        // the last argument is epsilon and is not used for polynomials
        _m2->positionDerivative(pMid, tmp->dt2dx, 1e-4);
        H.block(0, 0, _nPar1, 2) = tmp->h1 * tmp->dt2dx;
    } else
        _m1->transformPosAndErrors(where, pMid);
    if (_nPar2) {
        _m2->computeTransformAndDerivatives(pMid, outPoint, tmp->h2);
        H.block(_nPar1, 0, _nPar2, 2) = tmp->h2;
    } else
        _m2->transformPosAndErrors(pMid, outPoint);
}

/*! Sets the _nPar{1,2} and allocates H matrices accordingly, to
   avoid allocation at every call. If we did not care about dynamic
   allocation, we could just put the information of what moves and
   what doesn't into the SimpleGtransfoMapping. */
void TwoTransfoMapping::setWhatToFit(const bool fittingT1, const bool fittingT2) {
    if (fittingT1) {
        _nPar1 = _m1->getNpar();
        tmp->h1 = Eigen::MatrixX2d(_nPar1, 2);
    } else
        _nPar1 = 0;
    if (fittingT2) {
        _nPar2 = _m2->getNpar();
        tmp->h2 = Eigen::MatrixX2d(_nPar2, 2);
    } else
        _nPar2 = 0;
}

void TwoTransfoMapping::transformPosAndErrors(const FatPoint &where, FatPoint &outPoint) const {
    FatPoint pMid;
    _m1->transformPosAndErrors(where, pMid);
    _m2->transformPosAndErrors(pMid, outPoint);
}

void TwoTransfoMapping::positionDerivative(Point const &where, Eigen::Matrix2d &derivative,
                                           double epsilon) const {
    Eigen::Matrix2d d1, d2;  // seems that it does not trigger dynamic allocation
    _m1->positionDerivative(where, d1, 1e-4);
    FatPoint pMid;
    _m1->transformPosAndErrors(where, pMid);
    _m2->positionDerivative(pMid, d2, 1e-4);
    /* The following line is not a mistake. It is a consequence
       of chosing derivative(0,1) = d(y_out)/d x_in. */
    derivative = d1 * d2;
}

void TwoTransfoMapping::freezeErrorScales() {
    throw LSST_EXCEPT(pexExcept::TypeError,
                      " The routine  TwoTransfoMapping::freezeErrorScales() was thought to be useless and is "
                      "not implemented (yet)");
}
}  // namespace jointcal
}  // namespace lsst
