#include <iostream>
#include <fstream>
#include <iomanip>

#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/BaseStar.h"
#include "algorithm" // for copy

/* TO DO:
   think about imposing a maximum number of matches that may
   be discarded in Cleanup.
*/

namespace lsst {
namespace jointcal {



static double sq(double x) { return x*x;}

double StarMatch::Chi2(const Gtransfo &T) const
{
  FatPoint tr;
  T.TransformPosAndErrors(point1, tr);
  double vxx = tr.vx + point2.vx;
  double vyy = tr.vy + point2.vy;
  double vxy = tr.vxy + point2.vxy;
  double det = vxx*vyy-vxy*vxy;
  return (vyy*sq(tr.x-point2.x) + vxx*sq(tr.y-point2.y)
	  -2*vxy*(tr.x-point2.x)*(tr.y-point2.y))/det;
}


std::ostream& operator << (std::ostream &stream, const StarMatch &Match)
{
  stream << Match.point1.x << ' ' << Match.point1.y << ' '
	 << Match.point2.x << ' ' << Match.point2.y << ' '
	 << Match.distance << std::endl; return stream;
}

  std::ostream& operator << (std::ostream &stream, const StarMatchList &List)
{
  stream << " number of elements " << List.size() << std::endl;
  copy(List.begin(), List.end(), std::ostream_iterator<StarMatch>(stream));
  return stream;
}

static std::unique_ptr<double[]> chi2_array(const StarMatchList &L,const Gtransfo &T)
{
  unsigned s = L.size();
  auto res = std::unique_ptr<double[]>(new double[s]);
  unsigned count = 0;
  for (auto const &it: L)
    res[count++] = it.Chi2(T);
  return res;
}


static unsigned chi2_cleanup(StarMatchList &L,  const double Chi2Cut,
			     const Gtransfo &T)
{
  unsigned erased = L.RemoveAmbiguities(T);
  for (auto  smi = L.begin(); smi != L.end(); )
    {
      if (smi->chi2  > Chi2Cut)
	{
	  smi = L.erase(smi);
	  erased++;
	}
      else ++smi;
    }
  return erased;
}


#ifdef DO_WE_NEED_IT
static double *get_dist2_array(const StarMatchList &L, const Gtransfo &T)
{
  unsigned npair = L.size();
  if (npair == 0) return nullptr;
  std::unique_ptr<double[]> dist(new double [npair]);

  unsigned i=0;
  for (auto smi = L.begin(); smi != L.end(); ++smi, ++i)
    {
      const Point &p1 = smi->point1;
      const Point &p2 = smi->point2;
      dist[i] = p2.Dist2(T.apply(p1));
    }
  return dist;
}
#endif

int StarMatchList::Dof(const Gtransfo* T) const
{
  int npar = (T)? T->Npar() : (&(*transfo)?  transfo->Npar() : 0);
  int dof = 2*size() - npar;
  return (dof>0) ? dof : 0;
}




  /*! removes pairs beyond NSigmas in distance (where the sigma scale is
     set by the fit) and iterates until stabilization of the number of pairs.
     If the transfo is not assigned, it will be set to a GtransfoLinear. User
     can set an other type/degree using SetTransfo() before call. */
void StarMatchList::RefineTransfo(double NSigmas)
{
  double cut;
  unsigned nremoved;
  if (!transfo) transfo.reset(new GtransfoLin);
  do
    {
      int nused = size();
      if (nused <= 2) { chi2 = -1; break;}
      chi2 = transfo->fit(*this);
      /* convention of the fitted routines :
	 -  chi2 = 0 means zero degrees of freedom
              (this was not enforced in Gtransfo{Lin,Quad,Cub} ...)
	 -  chi2 = -1 means ndof <0 (and hence no possible fit)
	 --> in either case, refinement is over
         The fact that chi2 = 0 was not enforced when necessary means
	 that in this (rare) case, we were discarding matches at random....
         With GtransfoPoly::fit, this is no longer the case.
      */
      if (chi2 <= 0)  return;
      unsigned npair = int(size());
      if (npair == 0) break; // should never happen

      // compute some chi2 statistics
      std::unique_ptr<double[]> chi2_array(new double[npair]);
      unsigned count = 0;
      for (auto &starMatch: *this)
         chi2_array[count++] = starMatch.chi2 = starMatch.Chi2(*transfo);

      std::sort(chi2_array.get(), chi2_array.get()+npair);
      double median =  (npair&1)? chi2_array[npair/2] :
	(chi2_array[npair/2-1] + chi2_array[npair/2])*0.5;

      // discard outliers : the cut is understood as a "distance" cut
      cut = sq(NSigmas)*median;
      nremoved = chi2_cleanup(*this, cut, *transfo);
    }
  while (nremoved);
  dist2 = ComputeDist2(*this, *transfo);
}


/* not very robust : assumes that we went through Refine just before... */
double StarMatchList::Residual() const
{
  int deno = (2.*size() - transfo->Npar());
  return (deno>0) ? sqrt(dist2/deno) : -1; // is -1 a good idea?
}


void StarMatchList::SetDistance(const Gtransfo &Transfo)
{
  for (auto &smi: *this)
    smi.SetDistance(Transfo); // c'est compact
}


unsigned StarMatchList::RemoveAmbiguities(const Gtransfo &Transfo,
				     const int Which)
{
  if (!Which) return 0;
  SetDistance(Transfo);
  int initial_count = size();
  if (Which & 1)
    {
      sort(CompareS1);
      unique(SameS1);
    }
  if (Which & 2)
    {
      sort(CompareS2); unique(SameS2);
    }
  return (initial_count-size());
}


void StarMatchList::SetTransfoOrder(const int Order)
{
  if (Order==0) SetTransfo(std::make_shared<GtransfoLinShift>());
  else if (Order==1) SetTransfo(std::make_shared<GtransfoLin>());
  else SetTransfo(GtransfoPoly(Order));
  // might consider throwing if order does not make sense (e.g. >10)
  order = Order;
}



/* This routine should operate on a copy : RefineTransfo
   might shorten the list */
std::unique_ptr<Gtransfo> StarMatchList::InverseTransfo() /* it is not const although it tries not to change anything  */
{
  if (!transfo) return nullptr;

  auto old_transfo = transfo->Clone();
  double old_chi2 = chi2;

  Swap();
  SetTransfoOrder(order);
  RefineTransfo(3.);// keep same order
  auto inverted_transfo = transfo->Clone();
  SetTransfo(old_transfo.get());
  Swap();
  chi2 = old_chi2;

  return inverted_transfo;
}


void StarMatchList::CutTail(const int NKeep)
{
iterator si;
int count=0;
for (si = begin(); si != end() && count < NKeep; ++count, ++si);
erase(si, end());
}

void StarMatchList::Swap()
{
  for (auto &starMatch: *this)
    {
      starMatch.Swap() ;
    }
}

int StarMatchList::RecoveredNumber(double mindist) const
{
  int n = 0 ;
  GtransfoIdentity identity;
  for (auto const &starMatch: *this)
    {
      if (starMatch.Distance(identity) < mindist)
       n++ ;
    }
  return(n);
}


void StarMatchList::ApplyTransfo(StarMatchList &Transformed,
				 const Gtransfo *PriorTransfo,
				 const Gtransfo *PosteriorTransfo) const
{
  Transformed.clear();
  GtransfoIdentity id;
  const Gtransfo &T1 = (PriorTransfo)? *PriorTransfo : id;
  const Gtransfo &T2 = (PosteriorTransfo)? *PosteriorTransfo : id;

  for (auto const &starMatch: *this)
    {
      FatPoint p1;
      T1.TransformPosAndErrors(starMatch.point1,p1);
      FatPoint p2;
      T2.TransformPosAndErrors(starMatch.point2, p2);
      Transformed.push_back(StarMatch(p1, p2, starMatch.s1, starMatch.s2));
    }
}

void StarMatchList::SetChi2()
{
  chi2 = 0;
  for (auto &starMatch: *this)
    {
      starMatch.chi2 = starMatch.Chi2(*transfo);
      chi2 += starMatch.chi2;
    }
}

void StarMatchList::DumpTransfo(std::ostream &stream) const
{
  stream << " ================================================================" << std::endl
  << " Transformation between lists of order " << TransfoOrder() << std::endl
	 << *transfo //<< endl
  << " Chi2 = " << Chi2() << "  Residual = " << Residual() << std::endl
  << "  Number in the list = " << size() << std::endl
	 << " ================================================================" << std::endl;
}


double FitResidual(const double Dist2, const StarMatchList &S, const Gtransfo &T)
{
  return sqrt(Dist2/(2.*S.size()-T.Npar()));
}


double FitResidual(const StarMatchList &S, const Gtransfo &T)
{
  return FitResidual(ComputeDist2(S,T),S,T);
}

double ComputeDist2(const StarMatchList &S, const Gtransfo &T)
{
  double dist2 = 0;
  for (auto const &starMatch: S)
    dist2 += T.apply(starMatch.point1).Dist2(starMatch.point2);
  return dist2;
}

double computeChi2(const StarMatchList &L, const Gtransfo &T)
{
  unsigned s= L.size();
  std::unique_ptr<double[]> chi2s(chi2_array(L,T));
  double chi2 = 0;
  for (unsigned k=0; k<s; ++k) chi2 += chi2s[k];
  return chi2;
}


}} // end of namespaces
