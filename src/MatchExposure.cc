#include "lsst/jointcal/MatchExposure.h"
#include "lsst/jointcal/ExposureCatalog.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/ListMatch.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/Jointcal.h" // for JointcalControl
#include "lsst/pex/exceptions.h"

#include <fstream>
#include <vector>
#include <string>
using namespace std;

namespace lsst {
namespace jointcal {

/* If needed, documentation sits at the end of this file */


static Frame get_bounding_box(const ExposureStarList &L)
{
  // empty list ?
  if (L.begin() == L.end()) 
    {
      std::cout << " WARNING: get_bounding_box called with empty list "
		<< std::endl;
      return Frame();
    }
  const ExposureStar &first = **L.begin();
  double xmin= first.x;
  double ymin= first.y;
  double xmax = xmin;
  double ymax = ymin;
  for (auto i = L.begin(); i!= L.end(); ++i)
    {
      const Point &p = **i;
      if (p.x < xmin) xmin=p.x;
      if (p.x > xmax) xmax=p.x;
      if (p.y < ymin) ymin=p.y;
      if (p.y > ymax) ymax=p.y;
    }
  return Frame(xmin,ymin, xmax,ymax);
}

static bool accept_star_for_fit(const ExposureStar &S)
{
  return true;
}



// This is the routine that does the job
bool MatchExposure(ExposureCatalog &EC, const Point &TangentPoint, const JointcalControl &Control)
{
  // The arrangement has to be known by ExposureCatalog:
  const ChipArrangement &arrangement = EC.Arrangement(); 

  // assemble the catalog of the exposure (in tangent plane)
  ExposureStarList tpCat;
  EC.TangentPlaneCatalog(tpCat);

  // And get its boundaries, and take 10% margin.
  Frame coveredTpFrame = get_bounding_box(tpCat);
  Frame tpFrame= coveredTpFrame.Rescale(1.1);

  // the transformation that goes from the sky to the plane in which we match.
  // this is the tangent plane here, so:
  TanPix2RaDec tp2Sky(GtransfoLin(), TangentPoint);

  // OK, there are cases where this will not work, 
  // in particular close to the poles....
  Frame skyFrame = ApplyTransfo(tpFrame, tp2Sky, LargeFrame);

  // grab the reference catalog, from the provided file if any.
  BaseStarList refCat;
  UsnoCollect(skyFrame, tp2Sky, refCat);
  cout << "INFO: " << tpCat.size() << " objects in total in image" << endl;

  const std::vector<int> chips = EC.Chips();
  // Let us match now
  MatchConditions conditions;
	  
  conditions.NStarsL1 = 150;
  conditions.NStarsL2 = 150;
  conditions.MaxTrialCount = 10;
  conditions.SizeRatio = 1;
  // variations of plate scale due to focus are of the order of 1e-3.
  conditions.DeltaSizeRatio = 0.02; 
  cout << " INFO: running combinatorics: it takes a few seconds... " << endl;
  StarMatchList *match = MatchSearchRotShiftFlip((BaseStarList&)tpCat, refCat, conditions);
  if (!match) return false;

  GtransfoLin guessCorr = *dynamic_cast<const GtransfoLin *>(match->Transfo());
  cout << " INFO: correction to guessed WCS, found using " << match->size() 
       << " pairs" << endl<< guessCorr << endl;

  match = ListMatchCollect((BaseStarList&)tpCat, refCat, &guessCorr, Control.linMatchCut/3600);
  cout << " collected " << match->size() << " matches" << endl;

  int order = 2;
  match->SetTransfoOrder(order);
  match->RefineTransfo(3);

  GtransfoPoly guess(*dynamic_cast<const GtransfoPoly*>(match->Transfo()));
  // check if the whole FOV is within the window used to collect the 
  // reference data
  Frame actualTpFrame = ApplyTransfo(coveredTpFrame, guess,LargeFrame);
  double overlap = (actualTpFrame*tpFrame).Area()/actualTpFrame.Area();
  if (overlap<0.95)
    {
      skyFrame = ApplyTransfo(actualTpFrame, tp2Sky, LargeFrame);
      refCat.clear();
      cout << " INFO: too large a shift: we recollect the reference catalog " << endl;
      UsnoCollect(skyFrame, tp2Sky, refCat);
    }


  for (int iter = 0; iter<3; ++iter)
    {
      delete match;
      match = ListMatchCollect((BaseStarList&)tpCat, refCat, &guess, Control.secondMatchCut/3600);
      match->SetTransfoOrder(order);
      match->RefineTransfo(3);
      guess = *dynamic_cast<const GtransfoPoly*>(match->Transfo());
      cout << "re-collected " << match->size() << " matches (with refine)" << endl;
      cout << "new guess " << endl << guess << endl;
    }

  if (match->size() < Control.linMatchMinCount)
    {
      std::cout << "ERROR: not enough matches : giving up (if you think the match is OK, alter astromControl.linMatchMinCount"  << std::endl;
    
    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,"ERROR: not enough matches : giving up (if you think the match is OK, alter JointcalControl.linMatchMinCount");
    }

  // Now compute WCSs
  for (unsigned k=0; k<chips.size(); ++k)
    {
      unsigned chip = chips[k];
      int count = 0;
      StarMatchList subList;
      // extract from the large list and select on the fly
      for (StarMatchCIterator i = match->begin(); i!= match->end(); ++i)
	{
	  const ExposureStar &es = dynamic_cast<const ExposureStar &>(*(i->s1));
	  if (es.chip==chip)
	    {
	      count++;
	      if (accept_star_for_fit(es))
		subList.push_back(StarMatch(*es.original,*i->s2, es.original, &*(i->s2)));
	    }
	}
      cout << "INFO: before/after cleanup, number of matches for chip " << chip << " : " << count << '/' << subList.size() << endl; 

      if (subList.size() < Control.minMatchPerChip)
	{
	  std::cout << "ERROR: not enough matches for current chip: giving up (if you think the match is OK alter JointcalControl.minMatchPerChip)" 
		    << std::endl;
	  throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,"ERROR: not enough matches for current chip: giving up (if you think the match is OK alter JointcalControl.minMatchPerChip)");
	}
    

      double residual = 0;
      GtransfoPoly *pixToTP = NULL;// it can be set in 2 different ways;
      if (subList.size()>10) // we can carry out a fit
	{
	  unsigned match_order = Control.distortionDegree;
	  // ... if possible
	  while (((match_order+1)*(match_order+2))/2 > subList.size()) 
	    {match_order--;}
	  if (match_order != Control.distortionDegree) 
	    cout << "WARNING: we had to reduce distortion degree to " 
		 << match_order << " for chip  " << chip << endl;
	  pixToTP = new GtransfoPoly(match_order);

	  pixToTP->fit(subList);
	  residual = FitResidual(subList, *pixToTP)*3600; // in arcsec
	  cout << "INFO: 1d residual for distorted WCS :" << residual << " arcsec" << endl;
	}// we could actually fit something
      else // we could not fit: use the mapping
	{
	  const GtransfoPoly &corr = 
	    dynamic_cast<const GtransfoPoly &>(arrangement.Pix2TP(chip));
	  if (&corr==0) 
	    throw LSST_EXCEPT(lsst::pex::exceptions::TypeError,"ERROR: We cannot yet cook up a WCS if the \
arrangement file contains other stuff than polynomials. But there could be a way out.");
	  pixToTP = new GtransfoPoly(corr);
	}

#ifdef STORAGE
      // should we write?
      if (MatchPrefs.writeWCS == false) continue;
      ReducedImage ri(EC.ImageNames()[k]);
      string outFitsFileName(ri.FitsName());
      
      // where to output stuff : in the regular fits stuff or somewhere else? 
      if (MatchPrefs.wcsFileName != "") 
	{
	  if (strstr(MatchPrefs.wcsFileName.c_str(),"%s"))
	    {
	      char toto[256];
	      sprintf(toto, ("%s/"+MatchPrefs.wcsFileName).c_str(),
		      ri.Dir().c_str(), ri.Name().c_str());
	      outFitsFileName = toto;
	    }
	  else
	    outFitsFileName = ri.Dir()+"/"+MatchPrefs.wcsFileName;
	}
      if (outFitsFileName == ri.FitsName() && MatchPrefs.asciiWCS)
	{
	  cout << " ERROR : I refuse to overwrite a fits file by an ascii file : check the datacards " << std::endl;
	  cout << " no output " << endl;
	  return false;
	}
      
      cout << " writing WCS to " << outFitsFileName << endl; 
      {
	TanWCS2Header(outFitsFileName, *wcs);
	FitsHeader header(outFitsFileName, RW);
	header.AddOrModKey("RMATCHUS", residual, 
			   " 1d geom residual to ref catalog (arcsec)");
	header.AddOrModKey("NMATCHUS", int(subList.size()), 
			 " number of objects matched to ref catalog");
	string astromRef;
	if (MatchPrefs.astromCatalogName != "") 
	  astromRef = MatchPrefs.astromCatalogName;
	else 
	  {
	    char *usno_dir = getenv("USNODIR");
	    if (usno_dir) astromRef = BaseName(usno_dir); // strip path
	  }
	header.AddOrModKey("REFCAT",BaseName(astromRef), " Name of the ref. cat. astrom. and photom.");
      } // fitsfile is closed.

      if (MatchPrefs.asciiWCS) // convert to ASCII
	{
	  FitsHeader head(outFitsFileName);
	  unlink(outFitsFileName.c_str()); // means delete
	  ofstream s(outFitsFileName.c_str());
	  s << head << "END    " << std::endl;
	  s.close();
	}
#endif /* STORAGE */
    }// end loop on chips
  return true;
}

#ifdef STORAGE
static void usage(const char *pgname)
{
  cout << "usage : " << endl; 
  cout << pgname << ' ' << " <dbimage names> " << endl
       << " [-n] just print, won't alter WCS" << endl
       << " [-c <datacards> (example in datacards/match.datacards)]" << endl
       << " [-a <astrometric catalog>] (if applicable, superseeds value in datacards)" << endl;
  cout << " This code matches several CCDs from the same exposure " << endl
       << " to an astrometric catalogue (USNO-a by default) and" << endl
       << " sets WCSs (with distortions) in the fits headers" << endl;
  exit(1);
}


int main(int nargs, char** args)
{
  vector<string> images;
  try
    {
      for (int k=1; k<nargs; ++k) // first loop on arguments to read datacards name, if any
	{
	  const char *arg = args[k];
	  if (arg[0] != '-') continue;
	  if (arg[1] == 'c')
	    {
	      MatchPrefs.ReadCards(string(args[++k]));
	      break;
	    }
	}
      for (int k=1; k<nargs; ++k)
	{
	  const char *arg = args[k];
	  if (arg[0] == '-')
	    {
	      switch (arg[1])
		{
		case 'c' : k++; break; // datacards read above
		case 'n' : MatchPrefs.writeWCS = false; break;
		case 'a' : MatchPrefs.astromCatalogName = args[++k]; break;
		  if (access(MatchPrefs.astromCatalogName.c_str(),'r')!=0)
		    throw(PolokaException("ERROR: cannot find "+MatchPrefs.astromCatalogName));
		case 'h': 
		default : usage(args[0]);
		}
	    }
	  else
	    images.push_back(args[k]);
	} // end decoding arguments
      if (images.size() == 0) usage(args[0]);
      
      ExposureCatalog ec(images);
      if (match_exposure(ec)) return EXIT_SUCCESS;
    }
  catch(PolokaException p) {
    p.PrintMessage(cout);
  }    
  return EXIT_FAILURE;
}
#endif


/*
  This code matches several CCDs from the same exposure to a reference
  catalog.  In fact, it works poorly for a small numbers of CCDs,
  because the combinatorics is not done for that. This code is much
  faster than 36 invocations of matchusno, which should be used
  instead if only a few chips are to be matched. This code was
  initially developed to match very short exposures in order to
  identify without error the standards they contain. Matching one CCD at a
  time fails occasionally in these cases.

  How does it work?

  We rely on a set of exposures successfully matched from which we
  extract the mapping from pixels to the tangent plane. We average
  these mappings over a set of exposures, and store those. This is
  done by averagewcs.cc.

  These mappings are handled by ChipArrangement
  (chiparrangement.{cc,h}), which is a virtual class intended to be
  derived on a per-instrument basis.  It is in fact very simple code,
  and the scheme chosen to implement the mapping is not explicitely
  used in the matching code here. We have tested the code only on
  Megacam exposures, where a single ChipArrangement is adequate for
  matching exposures in all bands, over a wide range of airmasses,
  both before and after the lens flip. (In case you do not know, a
  lens of Megacam's wide field corrector was indeed flipped after
  ~1.5y of operation, with an improvement of IQ uniformity).

  Then, ExposureCatalog (exposurecatalog.{cc,h}) takes a set of
  DbImage's as input, figures out the instrument and grabs the
  corresponding ChipArrangement. It then transports the chip catalogs
  to the focal plane. The routine match_exposure (here) matches this
  ExposureCatalog to a reference catalog (USNO-a2) by default, and
  then fits separately the matches in each chip and writes the
  corresponding WCS's into the headers. Output can go into other files
  if required.

  The implementation of match_exposure (above) explicitely relies on
  the scheme chosen to implement the WCS distortions: it is the "PV"
  scheme used e.g. by swarp. If we want to make it more abstract, in
  order to accomodate other distortion schemes, this requires some
  changes in the Gtransfo class structures. We have to define a Wcs
  derived class which exposes the split between
  pixels-to-tangent-plane transformation and the de-projection
  transformation, and this class could be derived into different
  implementations. The trick we use to split between the CD matrix and
  the PV polynomials could be left to the concrete wcs classes, so
  that the code here becomes independent of the chosen WCS distortion.
  On the contrary, the code here does not rely on the chosen
  projection.  It turns out that at the moment, poloka only implements
  the RA___TAN stuff, so that this is purely rethoric.

  In the case some chips do not have enough matches to allow for a
  fit, we just assemble the pixel-to-tangent-plane mapping and the
  exposure match to cook up a WCS. This is indeed boot-strapping
  exposures to other exposures. At the moment, we rely on the input
  mappings being polynomials. We might consider cooking up a
  polynomial form the provided transformation if it is not the
  case. Not done yet.

  P.Astier 17/12/2014

 */


}} // end of namespaces
