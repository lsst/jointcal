#ifndef ASTROMFITCHI2__H
#define ASTROMFITCHI2__H

#include <string>
#include <iostream>
#include <sstream>

namespace lsst {
namespace meas {
namespace simastrom {

//! Simple structure to accumulate Chi2 and Ndof
struct Chi2
{
  double chi2;
  unsigned ndof;

  Chi2() : chi2(0), ndof(0) {};

  friend std::ostream& operator << (std::ostream& s, const Chi2 &C)
    {
      s << "Chi2/ndof : " << C.chi2 << '/' << C.ndof << '=' <<  C.chi2/C.ndof; return s;
    }


  //! this routine is the one called by the python print. 
  std::string __str__()
  {
    std::stringstream s;
    s << "Chi2/ndof : " << chi2 << '/' << ndof << '=' <<  chi2/ndof; 
    return s.str();
  }

  // Addentry has a third argument in order to make it compatible with an 
  //other stat accumulator.
  void AddEntry(double Inc, unsigned Dof, const void *M) 
  {chi2+= Inc; ndof += Dof;}

  void operator += (const Chi2 &R) 
  {chi2 += R.chi2; ndof += R.ndof;}

}; // end of struct Chi2




}}}
#endif /* ASTROMFIT__H */
