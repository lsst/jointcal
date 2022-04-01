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

#ifndef LSST_JOINTCAL_EIGENSTUFF_H
#define LSST_JOINTCAL_EIGENSTUFF_H

#include "lsst/pex/exceptions.h"

#include "Eigen/CholmodSupport"  // to switch to cholmod
#include "Eigen/Core"

using MatrixX2d = Eigen::Matrix<double, Eigen::Dynamic, 2>;

using SparseMatrixD = Eigen::SparseMatrix<double, 0, Eigen::Index>;

// To make our indices and triplets conform to Eigen's desire for taking a signed type
using IndexVector = std::vector<std::ptrdiff_t>;

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
    using Base = Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2>;
    using Base::m_cholmod;

public:
    using MatrixType = _MatrixType;
    using Index = typename MatrixType::Index;
    using RealScalar = typename MatrixType::RealScalar;

    CholmodSimplicialLDLT2() : Base() { init(); }

    CholmodSimplicialLDLT2(MatrixType const &matrix) : Base() {
        init();
        this->compute(matrix);
    }

    // this routine is the one we added
    void update(SparseMatrixD const &H, bool UpOrDown) {
        // check size
        Index const size = Base::m_cholmodFactor->n;
        EIGEN_UNUSED_VARIABLE(size);
        eigen_assert(size == H.rows());

        cholmod_sparse C_cs = viewAsCholmod(H);
        /* We have to apply the magic permutation to the update matrix,
        read page 117 of Cholmod UserGuide.pdf */
        // Using cholmod_l_* functions instead of cholmod_* because index is Eigen::Index instead of int.
        cholmod_sparse *C_cs_perm =
                cholmod_l_submatrix(&C_cs, (Eigen::Index*)Base::m_cholmodFactor->Perm,
                                  Base::m_cholmodFactor->n, nullptr, -1, true, true, &this->cholmod());
        assert(C_cs_perm);
        int isOk = cholmod_l_updown(UpOrDown, C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
        cholmod_l_free_sparse(&C_cs_perm, &this->cholmod());
        if (!isOk) {
            throw(LSST_EXCEPT(lsst::pex::exceptions::RuntimeError, "cholmod_update failed!"));
        }
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
