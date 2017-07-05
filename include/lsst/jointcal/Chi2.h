// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_CHI2_H
#define LSST_JOINTCAL_CHI2_H

#include <string>
#include <iostream>
#include <sstream>
#include <utility>

#include "lsst/jointcal/BaseStar.h"

static double sqr(double x) { return x * x; }

namespace lsst {
namespace jointcal {

/// Simple structure to accumulate chi2 and ndof
struct Chi2Statistic {
    double chi2;
    unsigned ndof;

    Chi2Statistic() : chi2(0), ndof(0){};

    friend std::ostream& operator<<(std::ostream& s, const Chi2Statistic& chi2) {
        s << "chi2/ndof : " << chi2.chi2 << '/' << chi2.ndof << '=' << chi2.chi2 / chi2.ndof;
        return s;
    }

    //! this routine is the one called by the python print.
    std::string __str__() {
        std::stringstream s;
        s << "Chi2/ndof : " << chi2 << '/' << ndof << '=' << chi2 / ndof;
        return s.str();
    }

    // Addentry has an ignored third argument in order to make it compatible with an
    // other stat accumulator.
    template <typename T>
    void addEntry(double inc, unsigned dof, T) {
        chi2 += inc;
        ndof += dof;
    }

    void operator+=(const Chi2Statistic& rhs) {
        chi2 += rhs.chi2;
        ndof += rhs.ndof;
    }

};  // end of struct Chi2

/*
 * A class to accumulate chi2 contributions together with pointers to the contributors.
 *
 * This structure lets one compute the chi2 statistics (average and variance) and directly point back
 * to the bad guys without relooping.
 * The Chi2Star routine makes it compatible with AstrometryFit's
 * accumulateStatImage and accumulateStatImageList.
 */
struct Chi2Star {
    double chi2;
    std::shared_ptr<BaseStar> star;

    Chi2Star(double chi2, std::shared_ptr<BaseStar> star) : chi2(chi2), star(std::move(star)) {}
    // for sorting
    bool operator<(const Chi2Star& rhs) const { return (chi2 < rhs.chi2); }
};

/// structure to accumulate chi2 and stars
struct Chi2List : public std::vector<Chi2Star> {
    void addEntry(double chi2, unsigned ndof, std::shared_ptr<BaseStar> star) {
        this->push_back(Chi2Star(chi2, std::move(star)));
    }

    /// Compute the average and std-deviation of these chisq values.
    std::pair<double, double> computeAverageAndSigma() {
        double sum = 0;
        double sum2 = 0;
        for (auto i : *this) {
            sum += i.chi2;
            sum2 += sqr(i.chi2);
        }
        double average = sum / this->size();
        double sigma = sqrt(sum2 / this->size() - sqr(average));
        return std::make_pair(average, sigma);
    }

    friend std::ostream& operator<<(std::ostream& s, const Chi2List& chi2List) {
        s << "chi2 per star : ";
        for (auto chi2 : chi2List) {
            s << *(chi2.star) << " chi2: " << chi2.chi2 << " ; ";
        }
        s << std::endl;
        return s;
    }
};

}  // namespace jointcal
}  // namespace lsst
#endif  // LSST_JOINTCAL_CHI2_H
