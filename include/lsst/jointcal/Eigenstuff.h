// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_EIGENSTUFF_H
#define LSST_JOINTCAL_EIGENSTUFF_H

#include "Eigen/CholmodSupport"  // to switch to cholmod

#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, 2> MatrixX2d;

typedef Eigen::SparseMatrix<double> SpMat;

/* Cholesky factorization class using cholmod, with the small-rank update capability.
 *
 * Class derived from Eigen's CholmodBase, to add the factorization
 * update capability to the interface. Besides this addition, it
 * behaves the same way as Eigen's native Cholesky factorization
 * classes. It relies on the simplicial LDLt factorization.
 *
 * @Seealso Eigen::CholmodSimplicialLDLT
 */
template <typename MatrixType, int UpLo = Eigen::Lower>
class CholmodSimplicialLDLT2
        : public Eigen::CholmodBase<MatrixType, UpLo, CholmodSimplicialLDLT2<MatrixType, UpLo>> {
    typedef Eigen::CholmodBase<MatrixType, UpLo, CholmodSimplicialLDLT2> Base;
    using Base::m_cholmod;

public:
    typedef MatrixType MatrixType;
    typedef typename MatrixType::Index Index;
    typedef typename MatrixType::RealScalar RealScalar;

    CholmodSimplicialLDLT2() : Base() { init(); }

    CholmodSimplicialLDLT2(MatrixType const &matrix) : Base() {
        init();
        this->compute(matrix);
    }

    // this routine is the one we added
    int update(SpMat const &h, bool upOrDown) {
        // check size
        Index const size = Base::m_cholmodFactor->n;
        EIGEN_UNUSED_VARIABLE(size);
        eigen_assert(size == h.rows());

        cholmod_sparse cCs = viewAsCholmod(h);
        /* We have to apply the magic permutation to the update matrix,
        read page 117 of Cholmod UserGuide.pdf */
        cholmod_sparse *cCsPerm =
                cholmod_submatrix(&cCs, (int *)Base::m_cholmodFactor->Perm, Base::m_cholmodFactor->n,
                                  nullptr, -1, true, true, &this->cholmod());
        assert(cCsPerm);
        int ret = cholmod_updown(upOrDown, cCsPerm, Base::m_cholmodFactor, &this->cholmod());
        cholmod_free_sparse(&cCsPerm, &this->cholmod());
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
