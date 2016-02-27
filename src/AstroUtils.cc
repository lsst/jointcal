#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <endian.h>
#include <string.h>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/Gtransfo.h"

namespace lsst {
namespace jointcal {

#ifndef M_PI
#define     M_PI            3.14159265358979323846  /* pi */
#endif

double RaStringToDeg(const std::string RaString)
{
int hours, minutes; double seconds;
if (sscanf(RaString.c_str(),"%d:%d:%lf", &hours, &minutes, &seconds) == 3)
  {
    return 15.*( double(hours) + double(minutes)/60. + seconds/3600.);
  }
// try to decode a decimal number
 char *endptr;
 const char* startptr = RaString.c_str();
 double ra = strtod(startptr, &endptr);
 if (endptr != startptr) // means successful conversion by strtod
   return ra;
 // should throw
 std::cerr << " cannot decode a right ascencion in " << RaString << std::endl;
 return 400.0; /* invalid value ! */
}

  // there is certainly a much more elegant way to do that.
  double DecStringToDeg(const std::string DecString)
{
  /* separators may be either : or ' */
  char dec[64];
  strcpy(dec, DecString.c_str());
  double minus_char = 1; /* no minus sign */
  for (char *p = dec; *p; p++) 
    {
      if (*p == ':' || *p == '\'' ) *p = ' ';
      if (*p == '-') minus_char = -1.;
    }
  int deg, minutes; double seconds;
  if (sscanf(dec,"%d %d %lf",&deg, &minutes, &seconds) == 3)
    {
      double sign = 1;
      if (deg < 0) sign = -1.;
      else if (deg == 0) sign = minus_char;
      return sign*(fabs(double(deg))+fabs(double(minutes))/60. + fabs(double(seconds))/3600.);
    }
  char *endptr;
  const char *startptr = DecString.c_str();
  double decVal = strtod(startptr, &endptr);
  if (startptr != endptr) // successful conversion bu strtod
    return decVal;
  // should throw
  std::cerr << " cannot decode a declination in : " << DecString << std::endl;
  return 200;
}


static double sq(const double &x) { return x*x;}



static inline void intswap(unsigned int &a_word)
{
  int tmp = be32toh(a_word);
  a_word = tmp;
}


static inline void read_a_star(FILE *ifp, unsigned int &raword, unsigned int  &decword, unsigned int &magword)
{ 
  fread(&raword,4,1,ifp);
  fread(&decword,4,1,ifp);
  fread(&magword,4,1,ifp);
  intswap(raword);
  intswap(decword);
  intswap(magword);
}





/****************************************************************
  Does the actual work of reading a USNO file. The format of the USNO
  files (picked up in the USNO distribution) is documented at the end
  of this source file
*****************************************************************/

//#define DEBUG_READ /* to actually debug the reading of a file */

static void readusno(const std::string &filebase,double minra,double maxra,
		     double mindec,double maxdec, UsnoColor Color,
		     BaseStarList &ApmList, 
		     const GtransfoLinShift& T= GtransfoLinShift(0,0))
{

  char line[80];
  int idx[96],len[96],i,raidx0,pos;
  float junk;
  FILE *ifp=NULL;
  unsigned int  raword,decword,magword;
  unsigned int minraword,maxraword,mindecword,maxdecword;
  double mag;
  
  std::string catname = std::string(filebase) + ".cat";
  std::string accname = std::string(filebase) + ".acc";



  /* Read the accelerator */
  
  if (!(ifp=fopen(accname.c_str(),"r"))) {
    std::cerr << "readusno error: Couldn't open " << accname << std::endl;
    return;
    }
  for (i=0 ; i<96 && !feof(ifp) ; ++i) {
    fgets(line,80,ifp);
    sscanf(line,"%f %d %d",&junk,&(idx[i]),&(len[i]));
    }
  fclose(ifp);
  if (i!=96) {
    std::cerr << "parsusno error: " << accname << " too short, only " << i << "  lines " << std::endl;
    return;
    }


/* Open the cat file and skip to the first index we'll need */

  if (!(ifp=fopen(catname.c_str(),"r"))) {
    std::cerr << "readusno error: Couldn't open " << catname << std::endl;
    return;
    }


raidx0=(int)floor(minra/3.75);       /* Indexed every 3.75deg = 15min */
pos=(idx[raidx0]-1)*12;
#ifdef DEBUG_READ
fprintf(stderr,"radix0=%d\n",raidx0);
fprintf(stderr,"Seeking to position %d in file %s\n",idx[raidx0],catname.c_str());
#endif
  if (fseek(ifp,pos,SEEK_SET) < 0) {
    std::cerr << "readusno: Error setting position in " << catname << std::endl;
    fclose(ifp);
    return ;
  }

/* Figure out the min and max ra and dec words */

  minraword=(unsigned int)floor(minra*3600*100 + 0.5);
  maxraword=(unsigned int)floor(maxra*3600*100 + 0.5);
  mindecword=(unsigned int)floor(mindec*3600*100 + 0.5);
  maxdecword=(unsigned int)floor(maxdec*3600*100 + 0.5);


  int mag_fact; /* mag = (magword/mag_fact)%1000, where mag_fact = 1 for R and 1000 for B  */
  if (Color == RColor) mag_fact = 1; 
  else if (Color == BColor) mag_fact = 1000;
  else 
    {
      std::cerr << " unknown color value requested in UsnoListGet " << std::endl; fclose(ifp); 
      return;
    }

#ifdef DEBUG_READ
fprintf(stderr,"minraword=%d, maxraword=%d\n",minraword,maxraword);
fprintf(stderr,"mindecword=%d, maxdecword=%d\n",mindecword,maxdecword);
#endif

/* Read until we find get into the RA range */

raword=0;

while (raword<minraword && !feof(ifp)) read_a_star(ifp, raword, decword, magword);

#ifdef DEBUG_READ
  fprintf(stderr,"Should be at minraword; minraword=%d, raword=%d\n",
  minraword,raword);
#endif

/* Read as long as we are in the right RA range; if a star read
   is in the right dec range, then keep it */

  while (raword<=maxraword  && !feof(ifp)) 
  {
    read_a_star(ifp, raword, decword, magword);
    if (decword>=mindecword && decword<=maxdecword && raword<=maxraword ) 
    {
      mag=(double)((magword/mag_fact)%1000)/10.;
      double ra = double(raword)/3600./100.;
      double dec = double(decword)/3600./100.-90;
      if (mag<99.9)
        {
          BaseStar *s = new BaseStar(ra,dec,mag);
	  /* We have to set errors here, because there is no single place 
	     below where it could easily be done. I wish we had a projection 
	     transfo in hand to define errors in the tangent plane.
	     For the value, we chose 0.3" r.m.s
	  */
	  s->vx = sq(0.3/3600/cos(dec*M_PI/180.));
	  s->vy = sq(0.3/3600);
	  s->vxy = 0;
	  if (s->x <minra || s->x > maxra) 
	    std::cout << " problem with usno star  : " << s << std::endl;
          *s = T.apply(*s); // no need to transfor errors : the transfo is a shift
	  ApmList.push_back(s);
        }        
    }
}
  
 fclose(ifp);
}


// read an ascii file containing alpha,delta mag, 1 item per line.
/* \page usno_file_format Astrometric file format
    You can provide your own match catalog to matchusno, either on the 
    command line, or via the USNOFILE environment variable. The code
    expects and ascii file and interprets the 3 first items of
    every line as alpha, delta (equatorial) and mag (where the mag 
    is mainly used to 
    isolate the brightest objects). Sexagesimal coordinates are accepted.
*/
  
static void read_ascii_astrom_file(const std::string &FileName, 
				  double MinRa,  double MaxRa,
				  double MinDec, double MaxDec, 
				  BaseStarList &ApmList, 
				  const GtransfoLinShift& T= GtransfoLinShift(0,0))
{
  FILE *file = fopen(FileName.c_str(),"r");
  if (!file)
    {
      std::cerr << " cannot open (supposedly ascii USNO file ) " 
		<< FileName << std::endl;
      return;
    }
  char line[512];
  while (fgets(line,512,file))
    {
      if (line[0] == '#') continue;
      char cra[32], cdec[32];
      double mag;
      if (sscanf(line,"%s %s %lf",cra,cdec,&mag) != 3)
	{
	  std::cerr << " cannot decode ra dec mag in :" << std::endl 
		    << line << std::endl;
	  break;
	}
      double ra = RaStringToDeg(cra);
      double dec = DecStringToDeg(cdec);
      if (( ra < MinRa) || (ra > MaxRa) || ( dec < MinDec) || ( dec > MaxDec))
	continue;
      BaseStar *s = new BaseStar(ra,dec,mag);
      *s = T.apply(*s);
      ApmList.push_back(s);
    }
  std::cout << " collected " << ApmList.size() 
	    << " objects from " << FileName << std::endl;
  fclose(file);
}


static void actual_usno_read(const std::string &usnodir,
			     double minra, double maxra, 
			     double mindec, double maxdec, UsnoColor Color,
			     BaseStarList &ApmList)
{
  int cat0,cat1;
  
  std::cout << " reading usno in window (" << minra << ',' 
	    << maxra << ") (" << mindec << ',' << maxdec << ")" << std::endl;

  /* Read parameters */
  
#ifdef DEBUG_READ
  fprintf(stderr,"Going to read USNO from directory %s\n",usnodir.c_str());
  fprintf(stderr,"minra=%lf, maxra=%lf\n",minra,maxra);
  fprintf(stderr,"mindec=%lf, maxdec=%lf\n",mindec,maxdec);
#endif

  /* Make sure the region is reasonable */

  if (minra == maxra || mindec == maxdec) return;
  if (minra > maxra) std::swap (minra,maxra);
  if (mindec > maxdec) std::swap(mindec,maxdec);
  if (mindec < -90.) mindec = -90;
  if (maxdec > 90.) maxdec = 90;

  /* Figure out which catalog files we need to read, and read them */

  mindec+=90.;
  maxdec+=90.;           /* Dec->spd */

  cat0=int(75*floor(mindec/7.5));
  cat1=int(75*floor(maxdec/7.5));
#ifdef DEBUG_READ
  fprintf(stderr,"cat0=%d,cat1=%d\n",cat0,cat1);
#endif
  for ( ; cat0<=cat1; cat0+=75) 
    {
      char cat_char[12];
      sprintf(cat_char,"zone%04d",cat0);
#ifdef DEBUG_READ
    fprintf(stderr,"Calling readusno on %s\n",cat_char);
#endif
    // some precaution if minra<0 and maxra >0 (have to read in 2 steps)
    std::string filename = usnodir+'/'+cat_char;
    if (minra<0 and maxra >0)
      {
	readusno(filename, minra+360.0,360.0,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(-360, 0));
	readusno(filename, 0,maxra,mindec,maxdec,Color, ApmList);
      }
    else if (minra<0 and maxra <0)
      {
	readusno(filename, minra+360.0,maxra+360.0,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(-360, 0));
      }
    else if (maxra>360. and minra<360.)
      {
	readusno(filename, minra, 360,mindec,maxdec,Color, ApmList);
	readusno(filename, 0,maxra-360,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(360,0));
      }
    else if (maxra>360. and minra>360.)
      {
	readusno(filename, minra-360.,maxra-360,mindec,maxdec,Color, ApmList, 
		 GtransfoLinShift(360,0));
      }
    else readusno(filename, minra,maxra,mindec,maxdec,Color, ApmList);
    }

#ifdef DEBUG_READ
  fprintf(stderr,"Returning.\n");
#endif
  int count = ApmList.size();
  std::cout << " collected " << count  << " objects" << std::endl;
}



//! ra's and dec's in degrees. Handles limits across alpha=0. (not with USNOFILE)
int UsnoRead(double minra, double maxra, 
	     double mindec, double maxdec, UsnoColor Color, 
	     BaseStarList &ApmList)
{
  const char *ascii_source = NULL;

  const char *env_var = getenv("USNOFILE");
  if (env_var) ascii_source = env_var;

  if (ascii_source)
    read_ascii_astrom_file(ascii_source, minra,maxra,mindec,maxdec, 
			   ApmList);
  else
    {
      const char *usno_dir = getenv("USNODIR");
      if (usno_dir)
	{
	  actual_usno_read(usno_dir, 
			   minra, maxra, mindec, maxdec, Color, 
			   ApmList);
	}
      else
	{
	  std::cerr << " ERROR : You should define USNODIR or USNOFILE env var, " 
		    << std::endl
		    << " or provide ASTROM_CATALOG_NAME via datacards to run this code " << std::endl;
	  return 0;
	}
    }
  return ApmList.size();
}


int UsnoRead(const Frame &W, UsnoColor Color, BaseStarList &ApmList)
{
  return UsnoRead(W.xMin, W.xMax, W.yMin, W.yMax, Color, ApmList);
}

Frame ApplyTransfo(const Frame& inputframe,const Gtransfo &T, const WhichTransformed W) 
{
  // 2 opposite corners
  double xtmin1, xtmax1, ytmin1, ytmax1;
  T.apply(inputframe.xMin,inputframe.yMin,xtmin1,ytmin1);
  T.apply(inputframe.xMax,inputframe.yMax,xtmax1,ytmax1);
  Frame fr1(std::min(xtmin1,xtmax1), std::min(ytmin1,ytmax1), 
	    std::max(xtmin1,xtmax1), std::max(ytmin1,ytmax1));
  // 2 other corners
  double xtmin2, xtmax2, ytmin2, ytmax2;
  T.apply(inputframe.xMin, inputframe.yMax, xtmin2, ytmax2);
  T.apply(inputframe.xMax, inputframe.yMin, xtmax2, ytmin2);
  Frame fr2(std::min(xtmin2,xtmax2), std::min(ytmin2,ytmax2), 
	    std::max(xtmin2,xtmax2), std::max(ytmin2,ytmax2));

  if (W == SmallFrame) return fr1*fr2;
  return fr1+fr2;
}



}} // end of namespaces
