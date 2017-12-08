#include <iostream>
#include <memory>
#include <cmath>   /* for floor */
#include <cstring> /* for memset*/

#include "lsst/log/Log.h"
#include "lsst/jointcal/Histo2d.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.Histo2d");
}

namespace lsst {
namespace jointcal {

Histo2d::Histo2d(int nnx, float mminx, float mmaxx, int nny, float mminy, float mmaxy) {
    nx = nnx;
    ny = nny;
    minx = mminx;
    miny = mminy;
    if (mmaxx != mminx)
        scalex = nx / (mmaxx - mminx);
    else {
        LOGL_WARN(_log, "Histo2d: minx = maxx requested");
        scalex = 1.0;
    }
    if (mmaxy != mminy)
        scaley = ny / (mmaxy - mminy);
    else {
        LOGL_WARN(_log, "Histo2d: maxy = miny requested");
        scaley = 1.0;
    }
    data = std::make_unique<float[]>(nx * ny);
    memset(data.get(), 0, nx * ny * sizeof(float));
}

Histo2d::Histo2d(const Histo2d &other) {
    memcpy(this, &other, sizeof(Histo2d));
    data = std::make_unique<float[]>(nx * ny);
    memcpy((data).get(), other.data.get(), nx * ny * sizeof(float));
}

bool Histo2d::indices(double x, double y, int &ix, int &iy) const {
    ix = (int)floor((x - minx) * scalex);
    if (ix < 0 || ix >= nx) return false;
    iy = (int)floor((y - miny) * scaley);
    return (iy >= 0 && iy < ny);
}

void Histo2d::fill(float x, float y, float weight) {
    int ix, iy;
    if (indices(x, y, ix, iy)) data[iy + ny * ix] += weight;
}

double Histo2d::maxBin(double &x, double &y) const {
    float *p, *pend;
    int imax = 0;
    float valmax = -1e30;

    for (p = data.get(), pend = p + nx * ny; pend - p; p++) {
        if (*p > valmax) {
            valmax = *p;
            imax = p - (data.get());
        }
    }
    int ix = imax / ny;
    int iy = imax - ix * ny;
    x = minx + ((float)ix + 0.5) / scalex;
    y = miny + ((float)iy + 0.5) / scaley;
    return valmax;
}

void Histo2d::zeroBin(double x, double y) {
    int ix, iy;
    if (indices(x, y, ix, iy)) data[iy + ny * ix] = 0;
}

double Histo2d::binContent(double x, double y) const {
    int ix, iy;
    if (indices(x, y, ix, iy)) return data[iy + ny * ix];
    LOGL_WARN(_log, "Histo2D::binContent outside limits requested");
    return -1e30;
}
}  // namespace jointcal
}  // namespace lsst
