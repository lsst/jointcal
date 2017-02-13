// This may look like C code, but it is really -*- C++ -*-

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/Frame.h"
namespace lsst {
namespace jointcal {


  class Frame;
  class Gtransfo;
  //! Transform a Frame through a Transfo.
  Frame ApplyTransfo(const Frame& inputframe,const Gtransfo &T,
		     const WhichTransformed W);

}} // end of namespaces
