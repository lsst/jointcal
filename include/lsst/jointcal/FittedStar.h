// -*- C++ -*-
#ifndef FITTEDSTAR__H
#define FITTEDSTAR__H

#include <iostream>
#include <fstream>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"

namespace lsst {
namespace jointcal {


class MeasuredStar;
class RefStar;
class Gtransfo;

/*! \file */


//! objects whose position is going to be fitted. Coordinates in Common Tangent Plane.

struct PmBlock
{
// proper motion in x and y. Units depend on how you fit them
  double pmx,pmy;
  double epmx, epmy, epmxy;
  double color; // OK it is unrelated, but associated in practice
  bool mightMove;


  PmBlock() : pmx(0), pmy(0), epmx(0), epmy(0), epmxy(0), color(0), mightMove(false) {};



};

//! The objects which have been measured several times. The MeasuredStar s measuring the same object in differenr CcdImage s point to the same FittedStar.
class FittedStar : public BaseStar, public PmBlock {

  friend class PhotomFit;
  friend class PhotomFit2;


  private:

  double mag;
  double emag;
  double col;
  int    gen;
  double wmag;
  unsigned indexInMatrix;
  int measurementCount;
  const RefStar *refStar;

  double flux2;
  double fluxErr;
  double fluxErr2;

 public:
  FittedStar() :
    BaseStar(), mag(-1), emag(-1), col(0.), gen(-1), wmag(0),
    indexInMatrix(-1), measurementCount(0), refStar(NULL),
    flux2(-1),
    fluxErr(-1),
    fluxErr2(-1) {}

  FittedStar(const BaseStar &B) :
    BaseStar(B), mag(-1), emag(-1), col(0.), gen(-1), wmag(0),
    indexInMatrix(0), measurementCount(0), refStar(NULL),
    flux2(-1),
    fluxErr(-1),
    fluxErr2(-1) {}

  //  FittedStar(const FittedStar& F)
  //    : BaseStar(F), mag(F.mag), emag(F.emag), col(F.col), gen(F.gen), wmag(F.wmag),
  //      index(F.index), measurementCount(F.measurementCount), refStar(F.refStar),
  //      flux2(F.flux2), fluxErr(F.fluxErr), fluxErr2(F.fluxErr2) {}

  //!
  FittedStar(const MeasuredStar &M);


  //!
  void ClearBeforeAssoc()
  {
    indexInMatrix = -1;
    measurementCount = 0;
    refStar = NULL;
    wmag = 0;
  }


  //!
  void dump(std::ostream & stream = std::cout) const
    { BaseStar::dumpn(stream);
    stream << " mcount "
	   << measurementCount << std::endl;
    }

  //!
  int MeasurementCount() const { return measurementCount;}

  //!
  int& MeasurementCount() { return measurementCount;}

  //! derived using available zero points in input images. In the absence ofZP, ZP= 0.
  double Mag() const { return mag;}
  double& Mag() { return mag; }
  double EMag() const { return emag;}
  double& EMag() { return emag; }
  double Col() const { return col;}
  double& Col() { return col; }
  int     Generation() const { return gen;}
  int&    Generation() { return gen; }

  //!
  void  SetMag(double Value) { mag = Value;}

  //! this routine will hopefully soon disappear.
  void AddMagMeasurement(const double MagValue,
			 const double MagWeight);

  //! index is a value that a fit can set and reread....
  void SetIndexInMatrix(const unsigned &Index){ indexInMatrix = Index;};

  //!
  int  IndexInMatrix() const { return indexInMatrix;}

  //!
  void SetRefStar(const RefStar*);

  //!
  const RefStar *GetRefStar() const { return refStar;};

  //! getters
  double         Flux() const { return flux; }
  double&        Flux() { return flux; }
  double         FluxErr() const { return fluxErr; }
  double&        FluxErr() { return fluxErr; }

  double         Flux2() const { return flux2; }
  double&        Flux2() { return flux2; }
  double         FluxErr2() const { return fluxErr2; }
  double&        FluxErr2() { return fluxErr2; }

  //! write stuff
  std::string       WriteHeader_(std::ostream& pr=std::cout, const char* i=NULL) const;
  virtual void      writen(std::ostream& s) const;
virtual void      read_it(std::istream& s, const char* format);
static BaseStar*  read(std::istream& s, const char* format);
};



class Gtransfo;

//! A list of FittedStar s. Such a list is typically constructed by Associations
class FittedStarList : public StarList<FittedStar>
{

  public :

  bool  inTangentPlaneCoordinates;

  //!
  FittedStarList() {inTangentPlaneCoordinates=true;}

#ifdef DO_WE_NEED_IT
  //! In the uncommon case case they are read from disk (written from a previous run)
    FittedStarList(const std::string &FileName);
#endif /* DO_WE_NEED_IT */

    /*! the routine that writes the output. TP2RaDec should convert Tp
      coordinates to sideral ones. GoodStars refer to those that are
      connected to ref Ccds and hence have a magnitude.
     */
    void WriteTuple(const std::string &FileName,
        const Gtransfo &TP2RaDec,
        const bool OnlyGoodStars = true);
};

typedef FittedStarList::const_iterator FittedStarCIterator;
typedef FittedStarList::iterator FittedStarIterator;
typedef CountedRef<FittedStar> FittedStarRef;



class FittedStarTuple {
  private :
  std::ofstream stream;

  public :
    FittedStarTuple( const std::string &FileName);

  //! out put some stuff for control plots.
  /*! F provides coordinates in tangent plane, RaDec the position on the sky */
  void AddEntry(const FittedStar &F, const Point &RaDec);

  ~FittedStarTuple() { stream.close();}
};

}} // end of namespaces

#endif /* FITTEDSTAR__H */
