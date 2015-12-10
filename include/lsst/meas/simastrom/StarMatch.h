// This may look like C code, but it is really -*- C++ -*-
#ifndef  STARMATCH__H 
#define  STARMATCH__H

#include <iostream>
#include <iterator> // for std::ostream_iterator
#include <algorithm> // for swap
#include <string>
#include <list>

#include "lsst/meas/simastrom/Point.h"
#include "lsst/meas/simastrom/CountedRef.h"
#include "lsst/meas/simastrom/BaseStar.h" // class definition used in inlined functions 
#include "lsst/meas/simastrom/Gtransfo.h" // inlined function calls Gtransfo::apply()

namespace lsst {
namespace meas {
namespace simastrom {



/*! \file
   \brief pairs of points

   Used to handles matches of star std::lists (image/image) or image/catalog
   to fit geometrical and photometric transformations.
   */

//! a pair of stars, usually belonging to different images. 
/*! One would normally assume that
      they are the same object of the sky. The object contains basically two 2d points (called
      later p1 and p2), and
      two pointers (unused by the class and its satellites), that enable in the end to trace back
      the stars in the caller data structures. */

//! A hanger for star associations
class StarMatch {
  
  friend class StarMatchList;

public: /* if one sets that private, then fitting routines will be a nightmare. we could set all the Transfo classes friend... */
  
  FatPoint point1, point2; //!< 2 points 
  CountedRef<const BaseStar> s1, s2; //!< the Star pointers (the pointer is in fact generic, pointed data is never used).
  double distance;
  double chi2;

public :
  //! constructor.
  /*! gives 2 points (that contain the geometry), plus pointers to the Star objects
    (which are there for user convenience). */
  StarMatch(const FatPoint &p1, const FatPoint &p2, const BaseStar *S1, const BaseStar *S2) : 
    point1(p1), point2(p2), s1(S1), s2(S2), distance(0.) {};
  
  // the next one would require that StarMatch knows BaseStar which is not mandatory for StarMatch to work
  // StarMatch(BaseStar *S1, BaseStar *S2) : point1(*S1), point2(*S2), s1(S1), s2(S2) {};

  //! returns the distance from T(p1) to p2. 
  double Distance(const Gtransfo &T) const { return  point2.Distance(T.apply(point1)); };

  //! returns the chi2 (using errors in the FatPoint's)
  double Chi2(const Gtransfo &T) const;

  //! to be used before sorting on distances. 
  void SetDistance(const Gtransfo &T) { distance = Distance(T);};
  //! returns the value computed by the above one.
  double Distance() const { return distance;}

  void Swap() { std::swap(point1, point2) ; std::swap(s1,s2) ; }

#ifndef SWIG
  friend bool DecFlux(const StarMatch & S1, const StarMatch & S2);
  friend bool IncreasingDistances(const StarMatch &one, const StarMatch &two);
  friend bool DecreasingDistances(const StarMatch &one, const StarMatch &two);
  friend bool DecPhoRatio(const StarMatch &S1, const StarMatch &S2);
#endif
  
  /* comparison that ensures that after a sort, duplicates are next one another */
  explicit StarMatch() {};
  
#ifndef SWIG
  friend std::ostream& operator << (std::ostream &stream, const StarMatch &Match); 
#endif

  ~StarMatch() {}


 private:
  //  bool operator <  (const StarMatch & other) const { return (s1 > other.s1) ? s1 > other.s1 : s2 > other.s2;}
 
  /* for unique to remove duplicates */
  bool operator == (const StarMatch &other) const { return (s1 == other.s1 && s2 == other.s2); };
  bool operator != (const StarMatch &other) const { return (s1 != other.s1 || s2 != other.s2); };

  friend bool CompareS1(const StarMatch &one, const StarMatch &two);
  friend bool SameS1(const StarMatch &one, const StarMatch &two);
  friend bool CompareS2(const StarMatch &one, const StarMatch &two);
  friend bool SameS2(const StarMatch &one, const StarMatch &two);
  
  //! enables \verbatim std::cout << mystarMatch << std::endl; \endverbatim

  //  ClassDef(StarMatch,1);
};


inline bool DecFlux(const StarMatch & S1, const StarMatch & S2)
{
  return(S1.s1->flux > S2.s1->flux); 
}

inline bool IncreasingDistances(const StarMatch &one, const StarMatch &two) 
{ 
  return(one.distance < two.distance);
}

inline bool DecreasingDistances(const StarMatch &one, const StarMatch &two) 
{
  return(one.distance > two.distance);
}

inline bool DecPhoRatio(const StarMatch &S1, const StarMatch &S2)
{
  return(S1.s1->flux/S1.s2->flux > S2.s1->flux/S2.s2->flux);
}


inline bool CompareS1(const StarMatch &one, const StarMatch &two)
{ 
  return ((one.s1 == two.s1) ? (one.distance < two.distance) : ( (const BaseStar*) one.s1 > (const BaseStar *) two.s1));
}

inline bool SameS1(const StarMatch &one, const StarMatch &two)
{ 
  return (one.s1 ==  two.s1);
}

inline bool CompareS2(const StarMatch &one, const StarMatch &two)
{ 
  return ((one.s2 == two.s2) ? (one.distance < two.distance) : ((const BaseStar *)one.s2 > (const BaseStar *) two.s2));
}

inline bool SameS2(const StarMatch &one, const StarMatch &two)
{ 
  return (one.s2 ==  two.s2);
}


/* =================================== StarMatchList ============================================================ */

#ifndef SWIG
std::ostream& operator << (std::ostream &stream, const StarMatch &Match);
#endif

//#ifdef TO_BE_FIXED
  typedef ::std::list<StarMatch>::iterator StarMatchIterator;
  typedef ::std::list<StarMatch>::const_iterator StarMatchCIterator;
//#endif

//! A std::list of star matches, 
/*! To be used as the argument to Gtransfo::fit routines. There is as
well a StarMatch::fit routine which fits a polynomial by default,
although the transfo may be user-provided. The
StarMatchList::RefineTransfo is a convenient tool to reject
outliers. Given two catalogs, one can assemble a StarMatchList using
utilities such as ListMatchCollect. StarMatchList's have write
capabilities.  NStarMatchList is a generalization of this 2-match to n-matches.
*/

#ifndef SWIG
std::ostream& operator << (std::ostream &stream, const StarMatchList &List);
#endif

class StarMatchList : public std::list<StarMatch> {

  private :
  int order;
  double chi2;
  double dist2;
  CountedRef<Gtransfo> transfo;



  public :

  void RefineTransfo(const double &NSigmas);


  //! enables to get a transformed StarMatchList. Only positions are transformed, not attached stars. const routine: "this" remains unchanged.
  void ApplyTransfo(StarMatchList &Transformed, 
		    const Gtransfo *PriorTransfo,
		    const Gtransfo *PosteriorTransfo = NULL) const;

  /* constructor */
  StarMatchList() : order(0), chi2(0){};

  //!carries out a fit with outlier rejection 


  //! enables to access the fitted transformation. 
  /*! Clone it if you want to store it permanently */
  const Gtransfo *Transfo() const { return &*(transfo);}

  //! access to the sum of squared residuals of the last call to RefineTransfo. 
  double Dist2() const { return dist2;}


  //! access to the chi2 of the last call to RefineTransfo. 
  double Chi2() const { return chi2;}


  //! returns the degree of freedom for the fit in x and y
  int Dof(const Gtransfo *T=NULL) const;

  //! returns the order of the used transfo 
  int TransfoOrder() const {return order;}

  
  //! swaps elements 1 and 2 of each starmatch in std::list. 
    void Swap () ;

  //! returns the average 1d Residual (last call to RefineTransfo)
    double Residual() const;

  /*! cleans up the std::list of pairs for pairs that share one of their stars, keeping the closest one.
     The distance is computed using Transfo. Which = 1 (2) removes ambiguities
     on the first (second) term of the match. Which=3 does both.*/
  unsigned RemoveAmbiguities(const Gtransfo &Transfo, const int Which=3);
  
  
  //! sets a transfo between the 2 std::lists and deletes the previous or default one.  No fit.
  void SetTransfo(const Gtransfo *Transfo) { transfo = Transfo->Clone();}
  //!
  void SetTransfo(const Gtransfo &Transfo) { transfo = Transfo.Clone();}

  //! set transfo according to the given order. 
  void SetTransfoOrder(const int Order);

  /*! returns the inverse transfo (Swap, fit(RefineTransfo) , and Swap). 
         The caller should delete the returned pointer. */
  Gtransfo *InverseTransfo();


  //! Sets the distance (residual) field of all std::list elements. Mandatory before sorting on distances 
  void SetDistance(const Gtransfo &Transfo);

  //! deletes the tail of the match std::list 
  void CutTail(const int NKeep);

  //! count the number of elements for which distance is < mindist
  int RecoveredNumber(double mindist) const;

  //! print the matching transformation quality (transfo, chi2, residual) 
  void DumpTransfo(std::ostream &stream = std::cout) const;

  ~StarMatchList() {/* should delete the transfo.... or use counted refs*/ };
  
  //!: without descriptor for l2tup 
  void write_wnoheader(std::ostream & pr=std::cout, const Gtransfo *T=NULL ) const ;


  //! write  StarMatchList with a header which  can be read by l2tup 
  void write(const std::string &filename, const Gtransfo *tf=NULL) const;
  void write(std::ostream &pr, const Gtransfo *tf=NULL) const;


    //! enables to dump a match std::list through : std::cout << List; 
  void dump(std::ostream & stream = std::cout) const { stream << *this; }

  //  ClassDef(StarMatchList,1)

    private:

  //! computes the chi2 even when there is no fit.
  void SetChi2();



  StarMatchList(const StarMatchList&); // copies nor properly handled 
  void operator=(const StarMatchList&);


};

//! r.m.s of 1 dim residual plots (corrected for fit d.o.f)
double FitResidual(const double Dist2, const StarMatchList &S, const Gtransfo &T);

//! r.m.s of 1 dim residual plots (corrected for fit d.o.f)
double FitResidual(const StarMatchList &S, const Gtransfo &T);

//! sum of distance squared			
double ComputeDist2(const StarMatchList &S, const Gtransfo &T);

//! the actual chi2
double ComputeChi2(const StarMatchList &L, const Gtransfo &T);

}}}
#endif /* STARMATCH__H */
