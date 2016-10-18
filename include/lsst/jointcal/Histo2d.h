#ifndef HISTO2D__H
#define HISTO2D__H


namespace lsst {
namespace jointcal {


class Histo2d {

 private:
  float *data;
  int nx,ny;
  float minx,miny;
  float scalex, scaley;

  bool indices(double X, double Y, int &ix, int &iy) const;

 public:

  Histo2d() {data=0;}
  Histo2d(int nx, float minx, float maxx, int ny, float miny,float maxy);

  Histo2d(const Histo2d &Other);

  void Fill(float x, float y, float weight=1.);

  double MaxBin(double &x, double &y) const ;

  void BinWidth(double &Hdx, double &Hdy) const { Hdx = 1./scalex; Hdy = 1./scaley;}

  double BinContent(double X, double Y) const;

  void ZeroBin(double X, double Y);

  ~Histo2d() { if (data) delete [] data;}

 private:
  void operator = (const Histo2d &Right);
};

}} // end of namespaces

#endif /* HISTO2D__H */
