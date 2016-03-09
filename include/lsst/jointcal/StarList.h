// This may look like C code, but it is really -*- C++ -*-
#ifndef STARLIST__H
#define STARLIST__H

#include <string>
#include <list>
#include <iostream>
#include <memory>

//#include "countedref.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/CountedRef.h"
//#include "lsst/jointcal/GlobalVal.h"

namespace lsst {
namespace jointcal {


class Frame;

//! std::lists of Stars. 
/*!     It is a template class, which means that the Star
type remains undefined until a user defines it. 
The std::list related operations (insertion,
sort, traversal) are to be carried out using STL 
list operations. Most of the
Star operations rely on routines to be provided in 
the Star class, usually
user defined. The instanciation of this class for 
BaseStar (i.e. the replacement 
of the formal parameter 'Star' by 'BaseStar') is 
called BaseStarList. 
Take care: what is stored is pointers on Star's and 
NOT Star's. This implies
that Stars being inserted in the std::list have to be 
obtained using 'new'.  */

/* CountedRef is very similar to std::shared_ptr, but with the latter
the following code does not work:
StarList<BaseStar> L;
for (auto i = L.begin(); i!= L.end(); ++i)
  {
    const BaseStar *s = *i;
    ...
  }
and the code is stuffed with such loops.
*/


  template<class Star> class StarList : public std::list <CountedRef<Star> > {
    //  GlobalVal glob;

public:
  typedef CountedRef<Star> Element;
  typedef typename std::list<Element>::const_iterator StarCIterator;
  typedef typename std::list<Element>::iterator StarIterator;


/* constructors */
//! : default constructor (empty std::list). 
  StarList() {};

 //! reads a StarList from a file, 
  /*!
     using the read method from the Star class. 
     See BaseStar for an example of implementation. */
  StarList(const std::string &FileName); 
  //! writes to a file
  /*!  calls iteratively the write method of the Star 
	class. It is unusable if the Star class does not 
	provide this functionnality.
	see BaseStar to see a possible implementation. 
  */


  Star *EmptyStar() const { return new Star();}

  int  write(const std::string &FileName) const;

  //! obvious meaning
  int read(const std::string &FileName);

#if (0)
  //! enables to access global values (lines starting with '@' in ascii files)
  GlobalVal &GlobVal() { return glob;}

  //! enables to access global values (lines starting with '@' in ascii files)
  const GlobalVal &GlobVal() const { return glob;}
#endif

/* destructor */
  virtual ~StarList() {};

  
  //! invokes dump(stream) for all Stars in the std::list. 
 void dump(std::ostream &stream = std::cout ) const { 
    for (auto p = this->begin(); 
	 p !=this->end(); ++p) (*p)->dump(stream);}

  //!a model routine to sort the std::list 
  /*! see DecreasingFlux() to see what it is, if you 
     want another sorting criterion) */
  // le premier de la std::liste a le plus grand flux
  void FluxSort();

  //! copy the head of the std::list at the  end of an other std::list (that may be empty on input)
  void ExtractHead(StarList<Star> &Out, int NHead) const;

  //! cuts the end of the std::list 
  void CutTail(const int NKeep);

  //! copy the part of the std::list which is included in the frame at the end of another std::list
  void ExtractInFrame(StarList<Star> &Out, const Frame &aFrame) const;
 //! cut the part of the std::list which is at a distance < mindist of the edges defined by frame.
  void CutEdges(const Frame &aFrame, float mindist);

  //! clears Copy and makes a copy of the std::list to Copy 
  void CopyTo(StarList<Star> &Copy) const;

  //! Clears the std::list 
  void ClearList() { CutTail(0);};

  //! enables to apply a geometrical transfo if Star is Basestar or derives from it.
  /*! could be extended to other type of transformations. */

  template<class Operator> void ApplyTransfo(const Operator &Op) 
  {for (auto p = this->begin(); p != this->end(); ++p) Op.TransformStar(*(*p));}


protected :
  int ascii_read(const std::string &FileName);

private :

  int read(std::istream & rd);


public :

  //!: with descriptor for l2tup 
  int  write(std::ostream & pr) const;


};

  //! enables \verbatim  std::cout << my_list; \endverbatim
#ifndef SWIG
template <class Star>  std::ostream & operator <<(std::ostream &stream, const StarList<Star> &List) 
    {List.dump(stream); return stream; }
#endif 

}} // end of namespaces

#endif /* STARLIST__H */


