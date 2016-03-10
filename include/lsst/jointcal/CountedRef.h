#ifndef COUNTEDREF__H
#define COUNTEDREF__H

//#include <cstdlib>
//#include <iostream>

#include "boost/intrusive_ptr.hpp"

namespace lsst {
namespace jointcal {


//! an implementation of "smart pointers" that counts references to an object. The object it "points" to has to derive from RefCount.

struct RefCount
{
  // the actual count has to be mutable so that one can pass a const reference.
  mutable unsigned refcount;


  RefCount() : refcount(0) {};
  RefCount(const RefCount &Other) : refcount (0) {};

};

/* The two follwing routines have to be provided by the end user when
   using boost::intrusive_ptr. So, here they are: */
template <class T> void intrusive_ptr_add_ref(T * p)
{
    if (p) p->refcount++;
}

template <class T> void intrusive_ptr_release(T* p)
{ // might check p->count>0
  if (p)
    if (--(p->refcount)==0) delete p;
}


// complicated typedef for templates (aka "alias typedef"):
template <typename T> using CountedRef = boost::intrusive_ptr<T> ;

}}

#endif /*COUNTEDREF__H */
