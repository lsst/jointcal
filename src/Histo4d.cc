#include <iostream>
#include <math.h> /* for floor */
#include <string.h> /* for memset*/
#include <algorithm> /* for sort */
#include <limits.h>

#include "lsst/jointcal/Histo4d.h"

namespace lsst {
namespace jointcal {

using namespace std;

SparseHisto4d::SparseHisto4d(const int N1, double Min1, double Max1,
			     const int N2, double Min2, double Max2,
			     const int N3, double Min3, double Max3,
			     const int N4, double Min4, double Max4,
			     const int nEntries)
{
  double indexMax = N1*N2*N3*N4;
  data = NULL;
  if (indexMax > double(INT_MAX))
    {
      cerr << " cannot hold a 4D histo with more than " << INT_MAX 	   << " values " <<  endl;
    }
  n[0] = N1;
  n[1] = N2;
  n[2] = N3;
  n[3] = N4;
  minVal[0] = Min1;
  minVal[1] = Min2;
  minVal[2] = Min3;
  minVal[3] = Min4;
  maxVal[0] = Max1;
  maxVal[1] = Max2;
  maxVal[2] = Max3;
  maxVal[3] = Max4;

  for (int i =0; i < 4; ++i)
    scale[i] = n[i]/(maxVal[i]-minVal[i]);
  data = new int[nEntries];
  dataSize = nEntries;
  ndata = 0;
  sorted = false;
}

int SparseHisto4d::code_value(const double X[4]) const
{
  int index = 0;
  for (int idim = 0; idim < 4 ; ++idim)
    {
      int  i = (int) floor(( X[idim] - minVal[idim])*scale[idim]);
      if (i <0 || i >= n[idim]) return - 1;
      index = index*n[idim] + i;
    }
  return index;
}

void SparseHisto4d::inverse_code(const int ACode, double X[4]) const
{
  int code = ACode;
  for (int i=3; i>=0; --i)
    {
      int bin = code%n[i];
      code /= n[i];
      X[i] = minVal[i] + ((double)bin + 0.5)/scale[i];
    }
}


void SparseHisto4d::sort()
{
  if (!sorted)
    {
      std::sort(data, data+ndata);
      sorted = true;
    }
}


void SparseHisto4d::Fill(const double X[4])

{
  int code = code_value(X);
  if (code <0) return;
  if (ndata == dataSize)
    {
      int* newData = new int[dataSize*2];
      memcpy(newData, data, dataSize*sizeof(data[0]));
      delete [] data;
      data = newData;
      dataSize *= 2;
    }
  data[ndata++] = code;
  sorted = false;
}

void SparseHisto4d::Fill(const double X1, const double X2, const double X3, const double X4)
{
  static double X[4];
  X[0] = X1; X[1] = X2; X[2] = X3; X[3] = X4; Fill(X);
}

int SparseHisto4d::MaxBin(double X[4])
{
  sort();
  if (ndata == 0) return 0;
  int maxval=data[0];
  int maxCount=1;
  int oldval = data[0];
  int currentCount = 1;
  for (int i=1; i < ndata; ++i)
    {
      if (data[i] == oldval) currentCount++;
      else currentCount = 1;
      if (currentCount > maxCount)
	{
	  maxCount = currentCount;
	  maxval = data[i];
	}
      oldval = data[i];
    }
  inverse_code(maxval, X);
  return maxCount;
}

void SparseHisto4d::ZeroBin(double X[4])
{
  sort();
  int code = code_value(X);
  // inefficient locator...
  int start = 0;
  while ((data[start] < code) && start < ndata) start++;
  // find how many identical entries we have
  int end = min(start+1,ndata);
  while (end < ndata && data[start] == data[end]) end++;
  int shift = end-start;
  int lastShift = ndata - (end-start);
  for (int i=start; i< lastShift; ++i) data[i] = data[i+shift];
  ndata -= shift;
}

void SparseHisto4d::BinLimits(const double X[4], const int Idim, double &Xmin, double &Xmax) const
{
  int code = code_value(X);
  double xCenter[4];
  inverse_code(code,xCenter);
  Xmin = xCenter[Idim] - 0.5/scale[Idim];
  Xmax = xCenter[Idim] + 0.5/scale[Idim];
}

void SparseHisto4d::dump() const
{
  for (int i=0; i<ndata; ++i) // DEBUG
    std::cout << data[i] << ' ';
  std::cout << std::endl;
}

}} // end of namespaces
