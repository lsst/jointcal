// This may look like C code, but it is really -*- C++ -*-
#ifndef USNOUTILS__H
#define USNOUTILS__H

#include "lsst/meas/simastrom/BaseStar.h"


namespace lsst {
namespace meas {
namespace simastrom {

enum UsnoColor { RColor , BColor};

  class Frame;

//!This routine reads the USNO catalog stars within the given bounds.
/*! The x and y coordinates of the given stars refer to RA and DEC respectively
expressed in  degrees. The location of the catalog should be set by the user
in the USNODIR environment variable. The catalog contains both
R and B magnitude. The Color argument has to be RColor or BColor. 
WARNING : The flux of the returned BaseStar's is in fact a magnitude. */


  int UsnoRead(const Frame &F, UsnoColor Color, BaseStarList &ApmList);

  //! accepts both sexagesimal and decimal values.
  double RaStringToDeg(const std::string RaString);

  double DecStringToDeg(const std::string DecString);


}}} // end of namespaces

#endif /* USNOUTILS__H */
