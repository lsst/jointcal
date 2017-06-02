// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_TRIPLETLIST_H
#define LSST_JOINTCAL_TRIPLETLIST_H

#include "Eigen/Sparse"

#include <vector>

namespace lsst {
namespace jointcal {

typedef Eigen::Triplet<double> Trip;

// at the moment this class implements the eigen format.
// it would be wise to implement it differently if talking to cholmod
class TripletList : public std::vector<Trip> {
private:
    unsigned _nextFreeIndex;

public:
    TripletList(int count) : _nextFreeIndex(0) { reserve(count); };

    void addTriplet(const unsigned i, const unsigned j, double val) { push_back(Trip(i, j, val)); }

    unsigned getNextFreeIndex() const { return _nextFreeIndex; }

    void setNextFreeIndex(unsigned index) { _nextFreeIndex = index; }
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_TRIPLETLIST_H
