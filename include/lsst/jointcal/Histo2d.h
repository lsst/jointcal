// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_HISTO2D_H
#define LSST_JOINTCAL_HISTO2D_H

namespace lsst {
namespace jointcal {

class Histo2d {
public:
    Histo2d() : _data() {}
    Histo2d(int nx, float minx, float maxx, int ny, float miny, float maxy);

    Histo2d(const Histo2d &other);

    void fill(float x, float y, float weight = 1.);

    double maxBin(double &x, double &y) const;

    void binWidth(double &hdx, double &hdy) const {
        hdx = 1. / scalex;
        hdy = 1. / scaley;
    }

    double binContent(double x, double y) const;

    void zeroBin(double x, double y);

    ~Histo2d() {}

private:
    void operator=(const Histo2d &right);
    bool _indices(double x, double y, int &ix, int &iy) const;

    std::unique_ptr<float[]> _data;
    int _nx, _ny;
    float _minx, _miny;
    float _scalex, _scaley;
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO2D_H
