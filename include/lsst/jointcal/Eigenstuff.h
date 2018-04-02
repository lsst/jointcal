// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_EIGENSTUFF_H
#define LSST_JOINTCAL_EIGENSTUFF_H

#include "Eigen/CholmodSupport"  // to switch to cholmod

#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, 2> MatrixX2d;

typedef Eigen::SparseMatrix<double> SparseMatrixD;

/* Cholesky factorization class using cholmod, with the small-rank update capability.
 *
 * Class derived from Eigen's CholmodBase, to add the factorization
 * update capability to the interface. Besides this addition, it
 * behaves the same way as Eigen's native Cholesky factorization
 * classes. It relies on the simplicial LDLt factorization.
 *
 * @Seealso Eigen::CholmodSimplicialLDLT
 */
template <typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSimplicialLDLT2
        : public Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2<_MatrixType, _UpLo>> {
    typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2> Base;
    using Base::m_cholmod;

public:
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;
    typedef typename MatrixType::RealScalar RealScalar;

    CholmodSimplicialLDLT2() : Base() { init(); }

    CholmodSimplicialLDLT2(MatrixType const &matrix) : Base() {
        init();
        this->compute(matrix);
    }

    // this routine is the one we added
    int update(SparseMatrixD const &H, bool UpOrDown) {
        // check size
        Index const size = Base::m_cholmodFactor->n;
        EIGEN_UNUSED_VARIABLE(size);
        eigen_assert(size == H.rows());

        cholmod_sparse C_cs = viewAsCholmod(H);
        /* We have to apply the magic permutation to the update matrix,
        read page 117 of Cholmod UserGuide.pdf */
        cholmod_sparse *C_cs_perm =
                cholmod_submatrix(&C_cs, (int *)Base::m_cholmodFactor->Perm, Base::m_cholmodFactor->n,
                                  nullptr, -1, true, true, &this->cholmod());
        assert(C_cs_perm);
        int ret = cholmod_updown(UpOrDown, C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
        cholmod_free_sparse(&C_cs_perm, &this->cholmod());
        assert(ret != 0);
        return ret;
    }

protected:
    void init() {
        m_cholmod.final_asis = 1;
        m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
        // In CholmodBase::CholmodBase(), the following statement is missing in
        // SuiteSparse 3.2.0.8. Fixed in 3.2.7
        Base::m_shiftOffset[0] = Base::m_shiftOffset[1] = RealScalar(0.0);
    }
};

#endif  // LSST_JOINTCAL_EIGENSTUFF_H
