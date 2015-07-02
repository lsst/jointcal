#ifndef HISTO4D__H
#define HISTO4D__H

namespace lsst {
namespace meas {
namespace simastrom {


//! A class to histogram in 4 dimensions. Uses Sparse storage. The number of bin is limited to 256 per dimension. Used in ListMatch.cc
class SparseHisto4d {

 private:
  int *data;
  int ndata; int dataSize;
  int n[4];
  double minVal[4],maxVal[4];
  double scale[4];
  bool sorted;
  
 public:
  SparseHisto4d() {}
  // obvious meanings. NEntries is used as the size of the primary allocation.
  SparseHisto4d(const int N1, double Min1, double Max1, 
		const int N2, double Min2, double Max2,
		const int N3, double Min3, double Max3,
		const int N4, double Min4, double Max4,
		const int NEntries);
  //!
  void Fill(const double X[4]);
  //!
  void Fill(const double X1, const double X2, const double X3, const double X4);
  //!
  int MaxBin(double X[4]);

  //!
  void ZeroBin(double X[4]);

  //! return the bin limits of dimension idim (0<=idim<4), around point X.
  void BinLimits(const double X[4], const int Idim, double &Xmin, double &Xmax) const;

  //!
  int NEntries() const { return ndata;}

  ~SparseHisto4d() { delete [] data;}

  // private:
  int code_value(const double X[4]) const;
  void inverse_code(const int ACode, double X[4]) const;
  void sort();
  void dump() const;

};

}}}

#endif /* HISTO2D__H */
