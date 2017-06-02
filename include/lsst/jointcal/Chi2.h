// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CHI2_H
#define LSST_JOINTCAL_CHI2_H

#include <string>
#include <iostream>
#include <sstream>

namespace lsst {
namespace jointcal {

//! Simple structure to accumulate Chi2 and Ndof
struct Chi2 {
    double chi2;
    unsigned ndof;

    Chi2() : chi2(0), ndof(0){};

    friend std::ostream& operator<<(std::ostream& s, const Chi2& chi2) {
        s << "chi2/ndof : " << chi2.chi2 << '/' << chi2.ndof << '=' << chi2.chi2 / chi2.ndof;
        return s;
    }

    //! this routine is the one called by the python print.
    std::string __str__() {
        std::stringstream s;
        s << "Chi2/ndof : " << chi2 << '/' << ndof << '=' << chi2 / ndof;
        return s.str();
    }

    // Addentry has a third argument in order to make it compatible with an
    // other stat accumulator.
    template <typename T>
    void addEntry(double inc, unsigned dof, T) {
        chi2 += inc;
        ndof += dof;
    }

    void operator+=(const Chi2& right) {
        chi2 += right.chi2;
        ndof += right.ndof;
    }

};  // end of struct Chi2
}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_CHI2_H
