// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_HISTO2D_H
#define LSST_JOINTCAL_HISTO2D_H

namespace lsst {
namespace jointcal {

class Histo2d {
public:
    Histo2d() : data() {}
    Histo2d(int nnx, float mminx, float mmaxx, int nny, float mminy, float mmaxy);

    Histo2d(Histo2d const &other);

    void fill(float x, float y, float weight = 1.);

    double maxBin(double &x, double &y) const;

    void binWidth(double &Hdx, double &Hdy) const {
        Hdx = 1. / scalex;
        Hdy = 1. / scaley;
    }

    double binContent(double x, double y) const;

    void zeroBin(double x, double y);

    ~Histo2d() {}

private:
    void operator=(Histo2d const &right);
    bool indices(double x, double y, int &ix, int &iy) const;

    std::unique_ptr<float[]> data;
    int nx, ny;
    float minx, miny;
    float scalex, scaley;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO2D_H
