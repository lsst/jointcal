#include "ndarray.h"
#include "Eigen/Core"

#include "lsst/afw/math/detail/TrapezoidalPacker.h"

#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

// ------------------ PhotometryTransfoChebyshev helpers ---------------------------------------------------

namespace {

// To evaluate a 1-d Chebyshev function without needing to have workspace, we use the
// Clenshaw algorith, which is like going through the recurrence relation in reverse.
// The CoeffGetter argument g is something that behaves like an array, providing access
// to the coefficients.
template <typename CoeffGetter>
double evaluateFunction1d(CoeffGetter g, double x, int size) {
    double b_kp2 = 0.0, b_kp1 = 0.0;
    for (int k = (size - 1); k > 0; --k) {
        double b_k = g[k] + 2 * x * b_kp1 - b_kp2;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    return g[0] + x * b_kp1 - b_kp2;
}

// This class imitates a 1-d array, by running evaluateFunction1d on a nested dimension;
// this lets us reuse the logic in evaluateFunction1d for both dimensions.  Essentially,
// we run evaluateFunction1d on a column of coefficients to evaluate T_i(x), then pass
// the result of that to evaluateFunction1d with the results as the "coefficients" associated
// with the T_j(y) functions.
struct RecursionArrayImitator {
    double operator[](int i) const {
        return evaluateFunction1d(coefficients[i], x, coefficients.getSize<1>());
    }

    RecursionArrayImitator(ndarray::Array<double const, 2, 2> const &coefficients_, double x_)
            : coefficients(coefficients_), x(x_) {}

    ndarray::Array<double const, 2, 2> coefficients;
    double x;
};

// Compute an affine transform that maps an arbitrary box to [-1,1]x[-1,1]
afw::geom::AffineTransform makeChebyshevRangeTransform(afw::geom::Box2D const &bbox) {
    return afw::geom::AffineTransform(
            afw::geom::LinearTransform::makeScaling(2.0 / bbox.getWidth(), 2.0 / bbox.getHeight()),
            afw::geom::Extent2D(-(2.0 * bbox.getCenterX()) / bbox.getWidth(),
                                -(2.0 * bbox.getCenterY()) / bbox.getHeight()));
}

// Initialize a "unit" Chebyshev
ndarray::Array<double, 2, 2> _identityChebyshev(size_t degree) {
    ndarray::Array<double, 2, 2> coeffs = ndarray::allocate(ndarray::makeVector(degree + 1, degree + 1));
    coeffs.deep() = 0.0;
    coeffs[0][0] = 1;
    return coeffs;
}
}  // namespace

PhotometryTransfoChebyshev::PhotometryTransfoChebyshev(size_t degree, afw::geom::Box2D const &bbox)
        : _bbox(bbox),
          _toChebyshevRange(makeChebyshevRangeTransform(bbox)),
          _coefficients(_identityChebyshev(degree)),
          _degree(degree),
          _nParameters((degree + 1) * (degree + 2) / 2) {}

PhotometryTransfoChebyshev::PhotometryTransfoChebyshev(ndarray::Array<double, 2, 2> const &coefficients,
                                                       afw::geom::Box2D const &bbox)
        : _bbox(bbox),
          _toChebyshevRange(makeChebyshevRangeTransform(bbox)),
          _coefficients(coefficients),
          _degree(coefficients.size() - 1),
          _nParameters((_degree + 1) * (_degree + 2) / 2) {}

double PhotometryTransfoChebyshev::transform(double x, double y, double instFlux) const {
    afw::geom::Point2D p = _toChebyshevRange(afw::geom::Point2D(x, y));
    return instFlux * evaluateFunction1d(RecursionArrayImitator(_coefficients, p.getX()), p.getY(),
                                         _coefficients.getSize<0>());
}

void PhotometryTransfoChebyshev::offsetParams(Eigen::VectorXd const &delta) {
    // NOTE: the indexing in this method and computeParameterDerivatives must be kept consistent!
    Eigen::VectorXd::Index k = 0;
    for (ndarray::Size j = 0; j <= _degree; ++j) {
        auto const iMax = _degree - j;  // to save re-computing `i+j <= degree` every inner step.
        for (ndarray::Size i = 0; i <= iMax; ++i, ++k) {
            _coefficients[j][i] -= delta[k];
        }
    }
}

namespace {
// The integral of T_n(x) over [-1,1]:
// https://en.wikipedia.org/wiki/Chebyshev_polynomials#Differentiation_and_integration
double integrateTn(int n) {
    if (n % 2 == 1)
        return 0;
    else
        return 2.0 / (1.0 - double(n * n));
}
}  // namespace

double PhotometryTransfoChebyshev::integrate() const {
    double result = 0;
    double determinant = _bbox.getArea() / 4.0;
    for (ndarray::Size j = 0; j < _coefficients.getSize<0>(); j++) {
        for (ndarray::Size i = 0; i < _coefficients.getSize<1>(); i++) {
            result += _coefficients[j][i] * integrateTn(i) * integrateTn(j);
        }
    }
    return result * determinant;
}

double PhotometryTransfoChebyshev::mean() const { return integrate() / _bbox.getArea(); }

void PhotometryTransfoChebyshev::computeParameterDerivatives(double x, double y, double instFlux,
                                                             Eigen::Ref<Eigen::VectorXd> derivatives) const {
    afw::geom::Point2D p = _toChebyshevRange(afw::geom::Point2D(x, y));
    // Algorithm: compute all the individual components recursively (since we'll need them anyway),
    // then combine them into the final answer vectors.
    Eigen::VectorXd Tnx(_degree + 1);
    Eigen::VectorXd Tmy(_degree + 1);
    Tnx[0] = 1;
    Tmy[0] = 1;
    if (_degree >= 1) {
        Tnx[1] = p.getX();
        Tmy[1] = p.getY();
    }
    for (ndarray::Size i = 2; i <= _degree; ++i) {
        Tnx[i] = 2 * p.getX() * Tnx[i - 1] - Tnx[i - 2];
        Tmy[i] = 2 * p.getY() * Tmy[i - 1] - Tmy[i - 2];
    }

    // NOTE: the indexing in this method and offsetParams must be kept consistent!
    Eigen::VectorXd::Index k = 0;
    for (ndarray::Size j = 0; j <= _degree; ++j) {
        auto const iMax = _degree - j;  // to save re-computing `i+j <= degree` every inner step.
        for (ndarray::Size i = 0; i <= iMax; ++i, ++k) {
            derivatives[k] = instFlux * Tmy[j] * Tnx[i];
        }
    }
}

Eigen::VectorXd PhotometryTransfoChebyshev::getParameters() const {
    Eigen::VectorXd parameters(_nParameters);
    // NOTE: the indexing in this method and offsetParams must be kept consistent!
    Eigen::VectorXd::Index k = 0;
    for (ndarray::Size j = 0; j <= _degree; ++j) {
        for (ndarray::Size i = 0; i + j <= _degree; ++i, ++k) {
            parameters[k] = _coefficients[j][i];
        }
    }

    return parameters;
}

}  // namespace jointcal
}  // namespace lsst
