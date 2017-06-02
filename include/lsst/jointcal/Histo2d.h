// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_HISTO2D_H
#define LSST_JOINTCAL_HISTO2D_H

namespace lsst {
namespace jointcal {

class Histo2d {
private:
    std::unique_ptr<float[]> data;
    int nx, ny;
    float minx, miny;
    float scalex, scaley;

    bool indices(double x, double y, int &ix, int &iy) const;

public:
    Histo2d() : data() {}
    Histo2d(int nx, float minx, float maxx, int ny, float miny, float maxy);

    Histo2d(const Histo2d &other);

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
    void operator=(const Histo2d &right);
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO2D_H
