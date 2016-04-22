#include <assert.h>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "lsst/jointcal/ChipArrangement.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/pex/exceptions.h"
namespace pexExcept = lsst::pex::exceptions;


using namespace std;

namespace lsst {
namespace jointcal {

//! expose this routine for documentation purposes
std::string ArrangementFileFromInstrumentName(const std::string &TelInstName)
{
  string ret(TelInstName);
  std::replace(ret.begin(), ret.end(),'/','_');
  return ret+"_average_wcs.txt";
}

PolyMappingArrangement::PolyMappingArrangement(const std::string &FileName, const unsigned NChips) : nchips(NChips)
{

  std::cout << "INFO : trying to read " << FileName << endl;
  
  auto t = ReadTransfoFile(FileName, tangentPlaneFrame);
  if (nchips && (t.size() != nchips))
    throw  LSST_EXCEPT(pex::exceptions::InvalidParameterError,
		       "ERROR: PolyMappingArrangement::PolyMappingArrangement\
 : could not read nchip mappings in "+FileName);
  
  if (!nchips) const_cast<unsigned &>(nchips) = t.size();

  pix2TP.resize(t.size());
  
  for (auto i = t.begin(); i != t.end(); ++i)
    {
      GtransfoPoly* x = dynamic_cast<GtransfoPoly*>(i->second.get()); // checking the types
      if (!x) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
				"ERROR: "+FileName+" does not contain GtransfoPoly mappings");
      pix2TP[i->first] = *x;
    }
}

const Gtransfo& PolyMappingArrangement::Pix2TP(const unsigned Chip) const
{
  if (Chip>=nchips)
    throw  LSST_EXCEPT(pex::exceptions::InvalidParameterError,
		       "PolyMappingArrangement::Pix2TP refers to a chip >= nchip");
  return pix2TP[Chip];
}


#ifdef STORAGE
const ChipArrangement* FindArrangement(const FitsHeader &Head)
{
  /* We are using TOADINST to find the file. If this turns out to be
     not practical (e.g. the mosaic evolved), we just have to add one
     more routine to VirtualInstrument (e.g. string
     ArrangementFileName() ) and only use the following scheme if the
     string comes back empty.
  */ 
  string telInstName = Head.KeyVal("TOADINST");
  unsigned nchips = NChips(Head);
  string fileName = ArrangementFileFromInstrumentName(telInstName);
  char *where = getenv("TOADSCARDS");
  if (where != 0) fileName = string(where)+"/"+fileName;
  else
    {
      cout << "FindArrangement searches files in $TOADSCARDS which is not defined" << endl;
      cout << "it usually points to poloka-core/datacards" << endl;
    }

  string errorMessage;
  ChipArrangement* ret=NULL;
  try {
    ret = new  PolyMappingArrangement(fileName, nchips);
  }
  catch (PolokaException &e) {errorMessage = e.message();};

  // When there are other types try them here. Error handling should be improved to distinguish ill-formed files from missing ones
  if (!ret)
    {
      std::string message(" FindArrangementFromName :"+telInstName+" "+errorMessage+". Have you run averagewcs ?");
      std::cout << message << std::endl;
      throw(PolokaException(message));
    }
  return ret;
}
#endif

		     


// identifier for the "frame line" in I/O's
#define FRAME_STRING "TPFrame"
#define CHIP_STRING "Chip"
 

// this routine just mirrors WriteTransfoFile. Change both accordingly
ChipTransfosType ReadTransfoFile(const string &FileName, Frame& TPFrame)
{
  ifstream s(FileName.c_str());
  if (!s) throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
			    " ReadTransfoFile : cannot open " + FileName);
  std::string frameString;
  s >>  frameString >> TPFrame.xMin >> TPFrame.yMin 
    >> TPFrame.xMax >> TPFrame.yMax;
  // sanity check
  if (frameString != FRAME_STRING || s.fail()) 
    throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "ReadTransfoFile : did not find the proper frame line ");
  // Passed. Read the transfos now
  unsigned count;
  s >> count;
  ChipTransfosType res;
  for (unsigned k=0; k<count; ++k)
    {
      std::string chipString;
      unsigned chip;
      s >> chipString >> chip;
      // sanity check
      if (chipString != CHIP_STRING || s.fail()) 
	throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "ReadTransfoFile : did not find the proper chip line ");

      res[chip] = GtransfoRead(s);
    }
   // Again a sanity check
   if (s.fail()) throw  LSST_EXCEPT(pex::exceptions::InvalidParameterError,"ReadTransfoFile : could not read all transfos");  
   s.close();
   return res;
}

/* Please, do not change the format here without changing the routine
   above.  consider that there are already such file around
   (check in ../datacards/ files like  *_average_wcs.txt) */

void WriteTransfoFile(const std::string &FileName, 
		      const ChipTransfosType &V, 
		      const Frame & TPFrame)
{
   ofstream s(FileName.c_str());
   if (!s) throw  LSST_EXCEPT(pex::exceptions::InvalidParameterError,
			      " write_transfos : cannot open " + FileName);
   s << FRAME_STRING 
     << ' ' << TPFrame.xMin << ' ' << TPFrame.yMin 
     << ' ' << TPFrame.xMax << ' ' << TPFrame.yMax << endl;
   s << V.size() << endl;
   for (auto i = V.begin(); i != V.end(); ++i)
     {
       s << CHIP_STRING << i->first << std::endl;
       i->second->Write(s);
     }
   // I am not sure if ofstream throws in case of serious problem.
   // but I am sure it raises a failbit:
   if (s.fail()) 
     throw  LSST_EXCEPT(pex::exceptions::InvalidParameterError,
			"write_transfos : could not write all transfos");  
   s.close(); // useless?
}




}} // end of namespaces 
