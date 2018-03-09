#include <iostream>
#include <cmath>     /* for floor */
#include <cstring>   /* for memset*/
#include <algorithm> /* for sort */
#include <memory>
#include <climits>

#include "lsst/log/Log.h"
#include "lsst/jointcal/Histo4d.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.Histo4d");
}

namespace lsst {
namespace jointcal {

SparseHisto4d::SparseHisto4d(const int n1, double min1, double max1, const int n2, double min2, double max2,
                             const int n3, double min3, double max3, const int n4, double min4, double max4,
                             const int nEntries) {
    double indexMax = n1 * n2 * n3 * n4;
    _data.reset();
    if (indexMax > double(INT_MAX)) {
        LOGLS_WARN(_log, "Cannot hold a 4D histo with more than " << INT_MAX << " values.");
    }
    _n[0] = n1;
    _n[1] = n2;
    _n[2] = n3;
    _n[3] = n4;
    _minVal[0] = min1;
    _minVal[1] = min2;
    _minVal[2] = min3;
    _minVal[3] = min4;
    _maxVal[0] = max1;
    _maxVal[1] = max2;
    _maxVal[2] = max3;
    _maxVal[3] = max4;

    for (int i = 0; i < 4; ++i) _scale[i] = _n[i] / (_maxVal[i] - _minVal[i]);
    _data = std::make_unique<int[]>(nEntries);
    _dataSize = nEntries;
    _ndata = 0;
    _sorted = false;
}

int SparseHisto4d::code_value(const double x[4]) const {
    int index = 0;
    for (int idim = 0; idim < 4; ++idim) {
        auto i = static_cast<int>(floor((x[idim] - _minVal[idim]) * _scale[idim]));
        if (i < 0 || i >= _n[idim]) return -1;
        index = index * _n[idim] + i;
    }
    return index;
}

void SparseHisto4d::inverse_code(int code, double x[4]) const {
    for (int i = 3; i >= 0; --i) {
        int bin = code % _n[i];
        code /= _n[i];
        x[i] = _minVal[i] + (static_cast<double>(bin) + 0.5) / _scale[i];
    }
}

void SparseHisto4d::sort() {
    if (!_sorted) {
        std::sort(_data.get(), _data.get() + _ndata);
        _sorted = true;
    }
}

void SparseHisto4d::fill(const double x[4])

{
    int code = code_value(x);
    if (code < 0) return;
    if (_ndata == _dataSize) {
        std::unique_ptr<int[]> newData(new int[_dataSize * 2]);
        memcpy(newData.get(), _data.get(), _dataSize * sizeof(_data[0]));
        _data.swap(newData);
        _dataSize *= 2;
    }
    _data[_ndata++] = code;
    _sorted = false;
}

void SparseHisto4d::fill(const double x1, const double x2, const double x3, const double x4) {
    static double x[4];
    x[0] = x1;
    x[1] = x2;
    x[2] = x3;
    x[3] = x4;
    fill(x);
}

int SparseHisto4d::maxBin(double x[4]) {
    sort();
    if (_ndata == 0) return 0;
    int maxval = _data[0];
    int maxCount = 1;
    int oldval = _data[0];
    int currentCount = 1;
    for (int i = 1; i < _ndata; ++i) {
        if (_data[i] == oldval) {
            currentCount++;
        } else {
            currentCount = 1;
        }
        if (currentCount > maxCount) {
            maxCount = currentCount;
            maxval = _data[i];
        }
        oldval = _data[i];
    }
    inverse_code(maxval, x);
    return maxCount;
}

void SparseHisto4d::zeroBin(double x[4]) {
    sort();
    int code = code_value(x);
    // inefficient locator...
    int start = 0;
    while ((_data[start] < code) && start < _ndata) start++;
    // find how many identical entries we have
    int end = std::min(start + 1, _ndata);
    while (end < _ndata && _data[start] == _data[end]) end++;
    int shift = end - start;
    int lastShift = _ndata - (end - start);
    for (int i = start; i < lastShift; ++i) _data[i] = _data[i + shift];
    _ndata -= shift;
}

void SparseHisto4d::binLimits(const double x[4], const int iDim, double &xMin, double &xMax) const {
    int code = code_value(x);
    double xCenter[4];
    inverse_code(code, xCenter);
    xMin = xCenter[iDim] - 0.5 / _scale[iDim];
    xMax = xCenter[iDim] + 0.5 / _scale[iDim];
}

void SparseHisto4d::dump() const {
    for (int i = 0; i < _ndata; ++i) {  // DEBUG
        std::cout << _data[i] << ' ';
    }
    std::cout << std::endl;
}
}  // namespace jointcal
}  // namespace lsst
