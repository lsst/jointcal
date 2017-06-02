#include <iostream>
#include <math.h>   /* for floor */
#include <string.h> /* for memset*/

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
    data.reset(new float[nx * ny]);
    memset(data.get(), 0, nx * ny * sizeof(float));
}

Histo2d::Histo2d(const Histo2d &Other) {
    memcpy(this, &Other, sizeof(Histo2d));
    data.reset(new float[nx * ny]);
    memcpy((this->data).get(), Other.data.get(), nx * ny * sizeof(float));
}

bool Histo2d::indices(double X, double Y, int &ix, int &iy) const {
    ix = (int)floor((X - minx) * scalex);
    if (ix < 0 || ix >= nx) return false;
    iy = (int)floor((Y - miny) * scaley);
    return (iy >= 0 && iy < ny);
}

void Histo2d::Fill(float X, float Y, float Weight) {
    int ix, iy;
    if (indices(X, Y, ix, iy)) data[iy + ny * ix] += Weight;
}

double Histo2d::MaxBin(double &X, double &Y) const {
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
    X = minx + ((float)ix + 0.5) / scalex;
    Y = miny + ((float)iy + 0.5) / scaley;
    return valmax;
}

void Histo2d::ZeroBin(double X, double Y) {
    int ix, iy;
    if (indices(X, Y, ix, iy)) data[iy + ny * ix] = 0;
}

double Histo2d::BinContent(double X, double Y) const {
    int ix, iy;
    if (indices(X, Y, ix, iy)) return data[iy + ny * ix];
    LOGL_WARN(_log, "Histo2D::BinContent outside limits requested");
    return -1e30;
}
}  // namespace jointcal
}  // namespace lsst
