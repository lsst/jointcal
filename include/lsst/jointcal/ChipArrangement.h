#ifndef CHIPARRANGEMENT__H
#define CHIPARRANGEMENT__H

#include <string>
#include <map>
//#include <vector>


#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/CountedRef.h"

namespace lsst {
namespace jointcal {


//! The class that exposes mappings from pixel coordinates to the tangent plane for a CCD mosaic. 
/*! It is meant to be derived for an actual mapping, such as the next classes. These mappings
are constructed from tools/averagewcs.cc, which can be consulted. */
class ChipArrangement
{
  std::string name;
  
 public :
  //!
  const std::string &InstrumentName() const {return name;}

  //! The actual transformations.
  virtual const Gtransfo& Pix2TP(const unsigned Chip) const = 0;
  //! Frame in the tangent plane w.r.t a conventional optical axis

  //! the size of the instrument.
  virtual Frame TangentPlaneFrame() const =0;


  virtual ~ChipArrangement() {};
};

//! A concrete implementation with polynomial mappings from CCD to TP coordinates
class PolyMappingArrangement : public ChipArrangement
{
  std::vector<GtransfoPoly> pix2TP;
  const unsigned nchips;
  Frame tangentPlaneFrame;
  
public:
  PolyMappingArrangement(const std::string &FileName, const unsigned NChips=0);

  virtual Frame TangentPlaneFrame() const
  {return tangentPlaneFrame;}
  
  const Gtransfo& Pix2TP(const unsigned Chip) const;
  
};




#ifndef SWIG

/* swig gets lost in namespaces for CountedRef*/
typedef std::map<unsigned,CountedRef<Gtransfo > > ChipTransfosType;

/* One way to describe the chip arrangement in the focal plane is to
map the Pixel coordinates of each chip on the tangent plane.  this is
described through transformations that can be written to and read from
files.  Here are the utilities for that: */

//! Read a set of transfos from a file. I guess it is useless to expose it
ChipTransfosType ReadTransfoFile(const std::string &FileName, Frame& TPFrame);

//! Write a set of transfos to a file.
void WriteTransfoFile(const std::string &FileName, const ChipTransfosType &V, const Frame & TPFrame); 

#endif /* SWIG */

}} // end of namespaces

#endif /* CHIPARRANGEMENT__H */

