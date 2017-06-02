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

    bool indices(double X, double Y, int &ix, int &iy) const;

public:
    Histo2d() : data() {}
    Histo2d(int nx, float minx, float maxx, int ny, float miny, float maxy);

    Histo2d(const Histo2d &Other);

    void Fill(float x, float y, float weight = 1.);

    double MaxBin(double &x, double &y) const;

    void BinWidth(double &Hdx, double &Hdy) const {
        Hdx = 1. / scalex;
        Hdy = 1. / scaley;
    }

    double BinContent(double X, double Y) const;

    void ZeroBin(double X, double Y);

    ~Histo2d() {}

private:
    void operator=(const Histo2d &Right);
};
}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_HISTO2D_H
