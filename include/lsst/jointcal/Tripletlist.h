#ifndef TRIPLETLIST__H
#define TRIPLETLIST__H

#include "Eigen/Sparse"

#include <vector>

namespace lsst {
namespace jointcal {


typedef Eigen::Triplet<double> Trip;

// at the moment this class implements the eigen format.
// it would be wise to implement it differently if talking to cholmod
class TripletList : public std::vector<Trip>
{
  unsigned nextFreeIndex;
 public :
  TripletList(int Count) {nextFreeIndex = 0; reserve(Count);};
    
  void AddTriplet(const unsigned i, const unsigned j, double val)
  {
    push_back(Trip(i,j,val));
  }

  unsigned NextFreeIndex() const
  {
    return nextFreeIndex;
  }

  void SetNextFreeIndex(unsigned Index)
  {
    nextFreeIndex = Index;
  }

};

}} // end of namespaces


#endif /* SPARSEALGEBRA__H */
