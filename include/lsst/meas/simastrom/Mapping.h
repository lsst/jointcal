#ifndef MAPPING__H
#define MAPPING__H


#include <vector>
#include "lsst/meas/simastrom/Eigenstuff.h"
#include "lsst/meas/simastrom/CountedRef.h"



namespace lsst {
namespace meas {
namespace simastrom {

class FatPoint;
class Point;

//! virtual class needed in the abstraction of the distortion model
class Mapping  : public RefCount
{

    public :
	
    //! Mumber of parameters in total
    virtual unsigned Npar() const = 0;


    //! Provides for this parameter set (of length Npar()) how they map into the "grand" fit
    virtual void GetMappingIndices(std::vector<unsigned> &Indices) const = 0;

    //! Actually applies the mapping and evaluates the derivatives w.r.t the fitted parameters.
    /*! This is grouped into a single call because for most models,
        evaluating the derivatives w.r.T parameters is not much longer
        than just transforming */
    virtual void ComputeTransformAndDerivatives(const FatPoint &Where,
						FatPoint &OutPos,
						Eigen::MatrixX2d &H) const = 0;
    //! The same as above but without the parameter derivatives (used to evaluate chi^2)
    virtual void TransformPosAndErrors(const FatPoint &Where,
				       FatPoint &OutPos) const = 0;

    //! Remember the error scale and freeze it
    //  virtual void FreezeErrorScales() = 0;

    virtual void OffsetParams(const double *Delta) = 0;

    //! The derivative w.r.t. position
    virtual void  PosDerivative(const Point &Where, Eigen::Matrix2d &Der, const double & Eps) const = 0;

    //!
    virtual ~Mapping() {};
    
};

// typedef std::list<CountedRef<Mapping> > MappingList;
// typedef MappingList::iterator MappingIterator;
// typedef MappingList::const_iterator MappingCIterator;

}}}
#endif /* MAPPING__H */ 
