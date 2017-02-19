#include <iostream>
#include <iomanip>
#include <algorithm>

#include "lsst/log/Log.h"
#include "lsst/jointcal/AstromFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/Mapping.h"

#include "lsst/jointcal/Gtransfo.h"
#include "Eigen/Sparse"
#include "Eigen/CholmodSupport" // to switch to cholmod
#include "lsst/pex/exceptions.h"
#include <fstream>
#include "lsst/jointcal/Tripletlist.h"

typedef Eigen::SparseMatrix<double> SpMat;


#if 0
template<typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSupernodalLLT2 : public Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSupernodalLLT2<_MatrixType, _UpLo> >
{
    typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSupernodalLLT2> Base;
    using Base::m_cholmod;

public:

    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;

    CholmodSupernodalLLT2() : Base() { init(); }
    //! The factorization happens in the constructor.
    CholmodSupernodalLLT2(const MatrixType& matrix) : Base()
    {
        init();
        this->compute(matrix);
    }

    //! this routine is the one we added.
    int update(const SpMat &H, const bool UpOrDown)
    {
        // check size
        const Index size = Base::m_cholmodFactor->n;
        EIGEN_UNUSED_VARIABLE(size);
        eigen_assert(size == H.rows());
        cholmod_sparse C_cs = viewAsCholmod(H);
        /* We have to apply the magic permutation to the update matrix,
           read page 117 of Cholmod UserGuide.pdf */
        cholmod_sparse *C_cs_perm = cholmod_submatrix(C_cs,
                                    Base::m_cholmodFactor->Perm,
                                    Base::m_cholmodFactor->n,
                                    nullptr, -1, nullptr, true, true,
                                    &this->cholmod());

        int ret = cholmod_updown(UpOrDown, &C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
        cholmod_free_sparse(C_cs_perm,  &this->cholmod());
        return ret;
    }

    ~CholmodSupernodalLLT2() {}
protected:
    void init()
    {
        m_cholmod.final_asis = 1;
        m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
        // In CholmodBase::CholmodBase(), the following statement is missing in
        // SuiteSparse 3.2.0.8. Fixed in 3.2.7
        Base::m_shiftOffset[0] = Base::m_shiftOffset[1] = RealScalar(0.0);
    }
};
#endif


//! Cholesky factorization class using cholmod, with the small-rank update capability.
/*! Class derived from Eigen's CholmodBase, to add the factorization
    update capability to the interface. Besides this addition, it
    behaves the same way as Eigen's native Cholesky factorization
    classes. It relies on the simplicial LDLt factorization.*/
template<typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSimplicialLDLT2 : public Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2<_MatrixType, _UpLo> >
{
    typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSimplicialLDLT2> Base;
    using Base::m_cholmod;

public:

    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Index Index;
    typedef typename MatrixType::RealScalar RealScalar;

    CholmodSimplicialLDLT2() : Base() { init(); }

    CholmodSimplicialLDLT2(const MatrixType& matrix) : Base()
    {
        init();
        this->compute(matrix);
    }

    // this routine is the one we added
    int update(const SpMat &H, const bool UpOrDown)
    {
        // check size
        const Index size = Base::m_cholmodFactor->n;
        EIGEN_UNUSED_VARIABLE(size);
        eigen_assert(size == H.rows());

        cholmod_sparse C_cs = viewAsCholmod(H);
        /* We have to apply the magic permutation to the update matrix,
        read page 117 of Cholmod UserGuide.pdf */
        cholmod_sparse *C_cs_perm = cholmod_submatrix(&C_cs,
                                    (int *) Base::m_cholmodFactor->Perm,
                                    Base::m_cholmodFactor->n,
                                    nullptr, -1, true, true,
                                    &this->cholmod());
        assert(C_cs_perm);
        int ret = cholmod_updown(UpOrDown, C_cs_perm, Base::m_cholmodFactor, &this->cholmod());
        cholmod_free_sparse(&C_cs_perm,  &this->cholmod());
        assert(ret != 0);
        return ret;
    }

    ~CholmodSimplicialLDLT2() {}
protected:
    void init()
    {
        m_cholmod.final_asis = 1;
        m_cholmod.supernodal = CHOLMOD_SIMPLICIAL;
        // In CholmodBase::CholmodBase(), the following statement is missing in
        // SuiteSparse 3.2.0.8. Fixed in 3.2.7
        Base::m_shiftOffset[0] = Base::m_shiftOffset[1] = RealScalar(0.0);
    }
};



using namespace std;

static double sqr(double x) {return x*x;}

//const double posErrorIncrement=0.02;

namespace {
    LOG_LOGGER _log = LOG_GET("jointcal.AstromFit");
}

namespace lsst {
namespace jointcal {

AstromFit::AstromFit(Associations &associations, DistortionModel *distortionModel, double posError) :
    _assoc(associations),  _distortionModel(distortionModel), _posError(posError)
{
    _LastNTrip = 0;
    _JDRef = 0;

    _posError = posError;

    _referenceColor = 0;
    _sigCol = 0;
    unsigned count = 0;
    for (auto const &i: _assoc.fittedStarList)
    {
        _referenceColor += i->color;
        _sigCol += sqr(i->color);
        count++;
    }
    if (count)
    {
        _referenceColor /= double(count);
        if (_sigCol > 0) _sigCol = sqrt(_sigCol/count - sqr(_referenceColor));
    }
    LOGLS_INFO(_log, "Reference Color: " << _referenceColor << " sig " << _sigCol);

    _nRefrac = _assoc.NBands();
    _refractionCoefficient = 0;

    _nMeasuredStars = 0;
    // The various _npar... are initialized in assignIndices.
    // Although there is no reason to adress them before one might be tempted by
    // evaluating a Chi2 rightaway, .. which uses these counts, so:
    assignIndices("");

}



#define NPAR_PM 2

/* ! this routine is used in 3 instances: when computing
the derivatives, when computing the Chi2, when filling a tuple.
*/
Point AstromFit::transformFittedStar(const FittedStar &fittedStar,
                                     const Gtransfo *sky2TP,
                                     const Point &refractionVector,
                                     double refractionCoeff,
                                     double mjd) const
{
    Point fittedStarInTP =  sky2TP->apply(fittedStar);
    if (fittedStar.mightMove)
    {
        fittedStarInTP.x += fittedStar.pmx*mjd;
        fittedStarInTP.y += fittedStar.pmy*mjd;
    }
    // account for atmospheric refraction: does nothing if color
    // have not been assigned
    // the color definition shouldbe the same when computing derivatives
    double color = fittedStar.color - _referenceColor;
    fittedStarInTP.x += refractionVector.x*color*refractionCoeff;
    fittedStarInTP.y += refractionVector.y*color*refractionCoeff;
    return fittedStarInTP;
}

/*! This is the first implementation of an error "model".  We'll
  certainly have to upgrade it. MeasuredStar provides the mag in case
  we need it.  */
static void tweakAstromMeasurementErrors(FatPoint &P, const MeasuredStar &Ms, double error)
{
    static bool called = false;
    static double increment = 0;
    if (!called)
    {
        increment = sqr(error); // was in Preferences
        called = true;
    }
    P.vx += increment;
    P.vy += increment;
}

// we could consider computing the chi2 here.
// (although it is not extremely useful)
void AstromFit::LSDerivatives1(const CcdImage &ccdImage,
                               TripletList &tList, Eigen::VectorXd &rhs,
                               const MeasuredStarList *msList) const
{
    /***************************************************************************/
    /**  Changes in this routine should be reflected into accumulateStatImage  */
    /***************************************************************************/
    /* Setup */
    /* this routine works in two different ways: either providing the
       ccdImage, of providing the MeasuredStarList. In the latter case, the
       ccdImage should match the one(s) in the list. */
    if (msList)
        assert ( (*(msList->begin()))->ccdImage == &ccdImage);

    // get the Mapping
    const Mapping *mapping = _distortionModel->getMapping(ccdImage);
    // count parameters
    unsigned npar_mapping = (_fittingDistortions) ? mapping->Npar() : 0;
    unsigned npar_pos = (_fittingPos) ? 2 : 0;
    unsigned npar_refrac = (_fittingRefrac) ? 1 : 0;
    unsigned npar_pm = (_fittingPM) ? NPAR_PM : 0;
    unsigned npar_tot =  npar_mapping + npar_pos + npar_refrac + npar_pm;
    // if (npar_tot == 0) this CcdImage does not contribute
    // any constraint to the fit, so :
    if (npar_tot == 0) return;
    vector<unsigned> indices(npar_tot, -1);
    if (_fittingDistortions)  mapping->GetMappingIndices(indices);

    // proper motion stuff
    double mjd = ccdImage.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = ccdImage.RefractionVector();
    // transformation from sky to TP
    const Gtransfo* sky2TP = _distortionModel->sky2TP(ccdImage);
    // reserve matrices once for all measurements
    GtransfoLin dypdy;
    // the shape of h (et al) is required this way in order to be able to
    // separate derivatives along x and y as vectors.
    Eigen::MatrixX2d h(npar_tot, 2), halpha(npar_tot, 2), hw(npar_tot, 2);
    Eigen::Matrix2d transW(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::VectorXd grad(npar_tot);
    // current position in the Jacobian
    unsigned kTriplets = tList.NextFreeIndex();
    const MeasuredStarList &catalog = (msList) ? *msList : ccdImage.getCatalogForFit();

    for (auto &i: catalog)
    {
        const MeasuredStar& ms = *i;
        if (!ms.IsValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = ms;
        tweakAstromMeasurementErrors(inPos, ms, _posError);
        h.setZero(); // we cannot be sure that all entries will be overwritten.
        FatPoint outPos;
        // should *not* fill h if whatToFit excludes mapping parameters.
        if (_fittingDistortions)
            mapping->ComputeTransformAndDerivatives(inPos, outPos, h);
        else mapping->TransformPosAndErrors(inPos, outPos);

        unsigned ipar = npar_mapping;
        double det = outPos.vx*outPos.vy - sqr(outPos.vxy);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, "Inconsistent measurement errors: drop measurement at " << Point(ms)
                       << " in image " << ccdImage.getName());
            continue;
        }
        transW(0, 0) = outPos.vy/det;
        transW(1, 1) = outPos.vx/det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy/det;
        // compute alpha, a triangular square root
        // of transW (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(transW(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = transW(0, 1)/alpha(0, 0);
        // DB - I think that the next line is equivalent to : alpha(1,1) = 1./sqrt(outPos.vy)
        // PA - seems correct !
        alpha(1, 1) = 1./sqrt(det*transW(0, 0));
        alpha(0, 1) = 0;

        const FittedStar *fs = ms.GetFittedStar();

        Point fittedStarInTP = transformFittedStar(*fs, sky2TP,
                               refractionVector,
                               _refractionCoefficient,
                               mjd);

        // compute derivative of TP position w.r.t sky position ....
        if (npar_pos > 0) // ... if actually fitting FittedStar position
        {
            sky2TP->Derivative(*fs, dypdy, 1e-3);
            // sign checked
            // TODO Still have to check with non trivial non-diagonal terms
            h(npar_mapping, 0) = -dypdy.A11();
            h(npar_mapping + 1, 0) = -dypdy.A12();
            h(npar_mapping, 1) = -dypdy.A21();
            h(npar_mapping + 1, 1) = -dypdy.A22();
            indices[npar_mapping] = fs->IndexInMatrix();
            indices.at(npar_mapping + 1) = fs->IndexInMatrix() + 1;
            ipar += npar_pos;
        }
        /* only consider proper motions of objects allowed to move,
        unless the fit is going to be degenerate */
        if (_fittingPM && fs->mightMove)
        {
            h(ipar, 0) = -mjd; // Sign unchecked but consistent with above
            h(ipar + 1, 1) = -mjd;
            indices[ipar] = fs->IndexInMatrix() + 2;
            indices[ipar + 1] = fs->IndexInMatrix() + 3;
            ipar += npar_pm;
        }
        if (_fittingRefrac)
        {
            /* if the definition of color changes, it has to remain
               consistent with transformFittedStar */
            double color = fs->color - _referenceColor;
            // sign checked
            h(ipar, 0) = -refractionVector.x*color;
            h(ipar, 1) = -refractionVector.y*color;
            indices[ipar] = _refracPosInMatrix;
            ipar += 1;
        }

        // We can now compute the residual
        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);

        // do not write grad = h*transW*res to avoid
        // dynamic allocation of a temporary
        halpha = h*alpha;
        hw = h*transW;
        grad = hw*res;
        // now feed in triplets and rhs
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar)
        {
            for (unsigned  ic = 0; ic < 2; ++ic)
            {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
                tList.AddTriplet(indices[ipar], kTriplets + ic, val);
#else
                tList.AddTriplet(kTriplets + ic, indices[ipar], val);
#endif
            }
            rhs(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2; // each measurement contributes 2 columns in the Jacobian
    } // end loop on measurements
    tList.SetNextFreeIndex(kTriplets);
}

// we could consider computing the chi2 here.
// (although it is not extremely useful)

#define HACK_REF_ERRORS 1. // used to isolate the measurement or ref terms

void AstromFit::LSDerivatives2(const FittedStarList &fsl, TripletList &tList, Eigen::VectorXd &rhs) const
{
    /* We compute here the derivatives of the terms involving fitted
       stars and reference stars. They only provide contributions if we
       are fitting positions: */
    if (! _fittingPos) return;
    /* the other case where the accumulation of derivatives stops
       here is when there are no RefStars */
    if (_assoc.refStarList.size() == 0) return;
    Eigen::Matrix2d w(2, 2);
    Eigen::Matrix2d alpha(2, 2);
    Eigen::Matrix2d h(2, 2), halpha(2, 2), hw(2, 2);
    GtransfoLin der;
    Eigen::Vector2d res, grad;
    unsigned indices[2 + NPAR_PM];
    unsigned kTriplets = tList.NextFreeIndex();
    /* We cannot use the spherical coordinates directly to evaluate
       Euclidean distances, we have to use a projector on some plane in
       order to express least squares. Not projecting could lead to a
       disaster around the poles or across alpha=0.  So we need a
       projector. We construct a projector and will change its
       projection point at every object */
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i: fsl)
    {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.SetTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.TransformPosAndErrors(*rs, rsProj);
        proj.Derivative(fs, der, 1e-4);
        // sign checked. TODO check that the off-diagonal terms are OK.
        h(0, 0) = -der.A11();
        h(1, 0) = -der.A12();
        h(0, 1) = -der.A21();
        h(1, 1) = -der.A22();
        // TO DO : account for proper motions.
        double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
        if (rsProj.vx <= 0 || rsProj.vy <= 0 || det <= 0)
        {
            LOGLS_WARN(_log, "RefStar error matrix not positive definite for:  " << *rs);
            continue;
        }
        w(0, 0) = rsProj.vy/det;
        w(0, 1) = w(1, 0) = -rsProj.vxy/det;
        w(1, 1) = rsProj.vx/det;
        w *= HACK_REF_ERRORS;
        det /= sqr(HACK_REF_ERRORS);
        // compute alpha, a triangular square root
        // of w (i.e. a Cholesky factor)
        alpha(0, 0) = sqrt(w(0, 0));
        // checked that  alpha*alphaT = transW
        alpha(1, 0) = w(0, 1)/alpha(0, 0);
        alpha(1, 1) = 1./sqrt(det*w(0, 0));
        alpha(0, 1) = 0;
        indices[0] = fs.IndexInMatrix();
        indices[1] = fs.IndexInMatrix() + 1;
        unsigned npar_tot = 2;
        /* TODO: account here for proper motions in the reference
        catalog. We can code the effect and set the value to 0. Most
        (all?)  catalogs do not even come with a reference epoch. Gaia
        will change that. When refraction enters into the game, one should
        pay attention to the orientation of the frame */

        /* The residual should be Proj(fs)-Proj(*rs) in order to be consistent
        with the measurement terms. Since P(fs) = 0, we have: */
        res[0] = -rsProj.x;
        res[1] = -rsProj.y;
        halpha = h*alpha;
        // grad = h*w*res
        hw = h*w;
        grad = hw*res;
        // now feed in triplets and rhs
        for (unsigned ipar = 0; ipar < npar_tot; ++ipar)
        {
            for (unsigned  ic = 0; ic < 2; ++ic)
            {
                double val = halpha(ipar, ic);
                if (val == 0) continue;
#if (TRIPLET_INTERNAL_COORD == COL)
                tList.AddTriplet(indices[ipar], kTriplets + ic, val);
#else
                tList.AddTriplet(kTriplets + ic, indices[ipar], val);
#endif
            }
            rhs(indices[ipar]) += grad(ipar);
        }
        kTriplets += 2; // each measurement contributes 2 columns in the Jacobian
    }
    tList.SetNextFreeIndex(kTriplets);
}


//! this routine computes the derivatives of all LS terms, including the ones that refer to references stars, if any
void AstromFit::LSDerivatives(TripletList &tList, Eigen::VectorXd &rhs) const
{
    auto L = _assoc.getCcdImageList();
    for (auto const &im: L)
    {
        LSDerivatives1(*im, tList, rhs);
    }
    LSDerivatives2(_assoc.fittedStarList, tList, rhs);
}


// This is almost a selection of lines of LSDerivatives1(CcdImage ...)
/* This routine (and the following one) is template because it is used
both with its first argument as "const CCdImage &" and "CcdImage &",
and I did not want to replicate it.  The constness of the iterators is
automagically set by declaring them as "auto" */

template <class ImType, class Accum>
void AstromFit::accumulateStatImage(ImType &image, Accum &accu) const
{
    /**********************************************************************/
    /**  Changes in this routine should be reflected into LSDerivatives1  */
    /**********************************************************************/
    /* Setup */
    // 1 : get the Mapping's
    const Mapping *mapping = _distortionModel->getMapping(image);
    // proper motion stuff
    double mjd = image.getMjd() - _JDRef;
    // refraction stuff
    Point refractionVector = image.RefractionVector();
    // transformation from sky to TP
    const Gtransfo* sky2TP = _distortionModel->sky2TP(image);
    // reserve matrix once for all measurements
    Eigen::Matrix2Xd transW(2, 2);

    auto &catalog = image.getCatalogForFit();
    for (auto const &i: catalog)
    {
        auto &ms = *i;
        if (!ms.IsValid()) continue;
        // tweak the measurement errors
        FatPoint inPos = ms;
        tweakAstromMeasurementErrors(inPos, ms, _posError);

        FatPoint outPos;
        // should *not* fill h if whatToFit excludes mapping parameters.
        mapping->TransformPosAndErrors(inPos, outPos);
        double det = outPos.vx*outPos.vy - sqr(outPos.vxy);
        if (det <= 0 || outPos.vx <= 0 || outPos.vy <= 0) {
            LOGLS_WARN(_log, " Inconsistent measurement errors :drop measurement at " << Point(ms)
                       << " in image " << image.getName());
            continue;
        }
        transW(0, 0) = outPos.vy/det;
        transW(1, 1) = outPos.vx/det;
        transW(0, 1) = transW(1, 0) = -outPos.vxy/det;

        const FittedStar *fs = ms.GetFittedStar();
        Point fittedStarInTP = transformFittedStar(*fs, sky2TP,
                               refractionVector,
                               _refractionCoefficient,
                               mjd);

        Eigen::Vector2d res(fittedStarInTP.x - outPos.x, fittedStarInTP.y - outPos.y);
        double chi2Val = res.transpose()*transW*res;

        accu.AddEntry(chi2Val, 2, &ms);
    }// end of loop on measurements
}


//! for a list of images.
template <class ListType, class Accum>
void AstromFit::accumulateStatImageList(ListType &list, Accum &accum) const
{
    for (auto &im: list)
    {
        accumulateStatImage(*im, accum);
    }
}

template <class Accum>
void AstromFit::accumulateStatRefStars(Accum &accum) const
{
    /* If you wonder why we project here, read comments in
       AstromFit::LSDerivatives2(TripletList &TList, Eigen::VectorXd &Rhs) */
    FittedStarList &fsl = _assoc.fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i: fsl)
    {
        FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.SetTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.TransformPosAndErrors(*rs, rsProj);
        // TO DO : account for proper motions.
        double rx = rsProj.x; // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
        double wxx = rsProj.vy/det;
        double wyy = rsProj.vx/det;
        double wxy = -rsProj.vxy/det;
        wxx *= HACK_REF_ERRORS;
        wyy *= HACK_REF_ERRORS;
        wxy *= HACK_REF_ERRORS;
        accum.AddEntry(wxx*sqr(rx) + 2*wxy*rx*ry + wyy*sqr(ry), 2, &fs);
    }
}

Chi2 AstromFit::computeChi2() const
{
    Chi2 chi2;
    accumulateStatImageList(_assoc.getCcdImageList(), chi2);
    // now ref stars:
    accumulateStatRefStars(chi2);
    // so far, ndof contains the number of squares.
    // So, subtract here the number of parameters.
    chi2.ndof -= _nParTot;
    return chi2;
}

//! a class to accumulate chi2 contributions together with pointers to the contributors.
/*! This structure allows to compute the chi2 statistics (average and
  variance) and directly point back to the bad guys without
  relooping. The Chi2Entry routine makes it compatible with
  accumulateStatImage and accumulateStatImageList. */
struct Chi2Entry
{
    double chi2;
    BaseStar *ps;

    Chi2Entry(double c, BaseStar *s): chi2(c), ps(s) {}
    // for sort
    bool operator < (const Chi2Entry &R) const {return (chi2 < R.chi2);}
};

struct Chi2Vect : public std::vector<Chi2Entry>
{
    void AddEntry(double Chi2Val, unsigned ndof, BaseStar *ps)
    { this->push_back(Chi2Entry(Chi2Val, ps));}

};

//! this routine is to be used only in the framework of outlier removal
/*! it fills the array of indices of parameters that a Measured star
    constrains. Not really all of them if you check. */
void AstromFit::getMeasuredStarIndices(const MeasuredStar &ms,
                                       std::vector<unsigned> &indices) const
{
    if (_fittingDistortions)
    {
        const Mapping *mapping = _distortionModel->getMapping(*ms.ccdImage);
        mapping->GetMappingIndices(indices);
    }
    const FittedStar *fs = ms.GetFittedStar();
    unsigned fsIndex = fs->IndexInMatrix();
    if (_fittingPos)
    {
        indices.push_back(fsIndex);
        indices.push_back(fsIndex + 1);
    }
    // For securing the outlier removal, the next block is just useless
    if (_fittingPM)
    {
        for (unsigned k = 0; k < NPAR_PM; ++k) indices.push_back(fsIndex + 2 + k);
    }
    /* Should not put the index of refaction stuff or we will not be
       able to remove more than 1 star at a time. */
}

//! contributions to derivatives of (presumambly) outlier terms. No discarding done.
void AstromFit::outliersContributions(MeasuredStarList &msOutliers,
                                      FittedStarList &fOutliers,
                                      TripletList &tList,
                                      Eigen::VectorXd &grad)
{
    // contributions from measurement terms:
    for (auto const &i: msOutliers)
    {
        MeasuredStar &out = *i;
        MeasuredStarList tmp;
        tmp.push_back(&out);
        const CcdImage &ccd = *(out.ccdImage);
        LSDerivatives1(ccd, tList, grad, &tmp);
    }
    LSDerivatives2(fOutliers, tList, grad);
}


//! Discards measurements and reference contributions contributing to the chi2 more than a cut, computed as <chi2>+nSigCut+rms(chi2) (statistics over contributions to the chi2). Returns the number of removed outliers. No refit done.
unsigned AstromFit::removeOutliers(double nSigCut,
                                   const std::string &measOrRef)
{
    MeasuredStarList msOutliers;
    FittedStarList fsOutliers;
    unsigned n = findOutliers(nSigCut, msOutliers, fsOutliers, measOrRef);
    removeMeasOutliers(msOutliers);
    removeRefOutliers(fsOutliers);
    return n;
}



//! Find Measurements and references contributing more than a cut, computed as <chi2>+nSigCut+rms(chi2). The outliers are NOT removed, and no refit is done.
/*! After returning from here, there are still measurements that
  contribute above the cut, but their contribution should be
  evaluated after a refit before discarding them. */
unsigned AstromFit::findOutliers(double nSigCut,
                                 MeasuredStarList &msOutliers,
                                 FittedStarList &fsOutliers,
                                 const std::string &measOrRef) const
{
    bool searchMeas = (measOrRef.find("Meas") != std::string::npos);
    bool searchRef = (measOrRef.find("Ref") != std::string::npos);

    // collect chi2 contributions
    Chi2Vect chi2s;
    chi2s.reserve(_nMeasuredStars + _assoc.refStarList.size());
    // contributions from measurement terms:
    if (searchMeas)
        accumulateStatImageList(_assoc.ccdImageList, chi2s);
    // and from reference terms
    if (searchRef)
        accumulateStatRefStars(chi2s);

    // do some stat
    unsigned nval = chi2s.size();
    if (nval == 0) return 0;
    sort(chi2s.begin(), chi2s.end());
    double median = (nval & 1) ? chi2s[nval/2].chi2 :
                    0.5*(chi2s[nval/2 - 1].chi2 + chi2s[nval/2].chi2);
    // some more stats. should go into the class if recycled anywhere else
    double sum = 0; double sum2 = 0;
    for (auto i = chi2s.begin(); i != chi2s.end(); ++i)
    {sum += i->chi2; sum2 += sqr(i->chi2);}
    double average = sum/nval;
    double sigma = sqrt(sum2/nval - sqr(average));
    LOGLS_DEBUG(_log, "RemoveOutliers chi2 stat: mean/median/sigma "
               << average << '/' << median << '/' << sigma);
    double cut = average + nSigCut*sigma;
    /* For each of the parameters, we will not remove more than 1
       measurement that contributes to constraining it. Keep track using
       of what we are touching using an integer vector. This is the
       trick that Marc Betoule came up to for outlier removals in "star
       flats" fits. */
    Eigen::VectorXi affectedParams(_nParTot);
    affectedParams.setZero();

    unsigned nOutliers = 0; // returned to the caller
    // start from the strongest outliers.
    for (auto i = chi2s.rbegin(); i != chi2s.rend(); ++i)
    {
        if (i->chi2 < cut) break; // because the array is sorted.
        vector<unsigned> indices;
        indices.reserve(100); // just there to limit reallocations.
        /* now, we want to get the indices of the parameters this chi2
        term depends on. We have to figure out which kind of term it
         is; we use for that the type of the star attached to
         the Chi2Entry. */
        MeasuredStar *ms = dynamic_cast<MeasuredStar *>(i->ps);
        FittedStar *fs = nullptr;
        if (!ms) // it is reference term.
        {
            fs = dynamic_cast<FittedStar *>(i->ps);
            indices.push_back(fs->IndexInMatrix());
            indices.push_back(fs->IndexInMatrix() + 1); // probably useless
            /* One might think it would be useful to account for PM
               parameters here, but it is just useless */
        }
        else // it is a measurement term.
        {
            getMeasuredStarIndices(*ms, indices);
        }

        /* Find out if we already discarded a stronger outlier
        constraining some parameter this one constrains as well. If
         yes, we keep this one, because this stronger outlier could be
         causing the large chi2 we have in hand.  */
        bool drop_it = true;
        for (auto const &i: indices)
            if (affectedParams(i) != 0) drop_it = false;

        if (drop_it) // store the outlier in one of the lists:
        {
            if (ms) // measurement term
                msOutliers.push_back(ms);
            else // ref term
                fsOutliers.push_back(fs);
            // mark the parameters as directly changed when we discard this chi2 term.
            for (auto const &i: indices)
                affectedParams(i)++;
            nOutliers++;
        }
    } // end loop on measurements/references
    LOGLS_INFO(_log, "findOutliers: found "
               << msOutliers.size() << " meas outliers and "
               << fsOutliers.size () << " ref outliers ");

    return nOutliers;
}


void AstromFit::removeMeasOutliers(MeasuredStarList &outliers)
{
    for (auto &i: outliers)
    {
        MeasuredStar &ms = *i;
        FittedStar *fs = const_cast<FittedStar *>(ms.GetFittedStar());
        ms.SetValid(false);
        fs->MeasurementCount()--; // could be put in SetValid
    }
}


void AstromFit::removeRefOutliers(FittedStarList &outliers)
{
    for (auto &i: outliers)
    {
        FittedStar &fs = *i;
        fs.setRefStar(nullptr);
    }
}

void AstromFit::assignIndices(const std::string &whatToFit)
{
    _whatToFit = whatToFit;
    LOGLS_INFO(_log, "assignIndices: Now fitting " << whatToFit);
    _fittingDistortions = (_whatToFit.find("Distortions") != string::npos);
    _fittingPos = (_whatToFit.find("Positions") != string::npos);
    _fittingRefrac = (_whatToFit.find("Refrac") != string::npos);
    if (_sigCol == 0 && _fittingRefrac)
    {
        LOGLS_WARN(_log, "Cannot fit refraction coefficients without a color lever arm. Ignoring refraction.");
        _fittingRefrac = false;
    }
    _fittingPM = (_whatToFit.find("PM") != string::npos);
// When entering here, we assume that whatToFit has already been interpreted.


    _nParDistortions = 0;
    if (_fittingDistortions)
        _nParDistortions = _distortionModel->assignIndices(0, _whatToFit);
    unsigned ipar = _nParDistortions;

    if (_fittingPos)
    {
        FittedStarList &fsl = _assoc.fittedStarList;
        for (auto const &i: fsl)
        {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when updating (offsetParams())
            // - in GetMeasuredStarIndices
            fs.SetIndexInMatrix(ipar);
            ipar += 2;
            if ((_fittingPM) & fs.mightMove) ipar += NPAR_PM;
        }
    }
    unsigned _nParPositions = ipar - _nParDistortions;
    if (_fittingRefrac)
    {
        _refracPosInMatrix = ipar;
        ipar += _nRefrac;
    }
    _nParTot = ipar;
}



void AstromFit::offsetParams(const Eigen::VectorXd& delta)
{
    if (delta.size() != _nParTot)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "AstromFit::offsetParams : the provided vector length is not compatible with the current whatToFit setting");
    if (_fittingDistortions)
        _distortionModel->offsetParams(delta);

    if (_fittingPos)
    {
        FittedStarList &fsl = _assoc.fittedStarList;
        for (auto const &i: fsl)
        {
            FittedStar &fs = *i;
            // the parameter layout here is used also
            // - when filling the derivatives
            // - when assigning indices (assignIndices())
            unsigned index = fs.IndexInMatrix();
            fs.x += delta(index);
            fs.y += delta(index + 1);
            if ((_fittingPM) & fs.mightMove)
            {
                fs.pmx += delta(index + 2);
                fs.pmy += delta(index + 3);
            }
        }
    }
    if (_fittingRefrac)
    {
        _refractionCoefficient += delta(_refracPosInMatrix);
    }
}

// should not be too large !
#ifdef STORAGE
static void write_sparse_matrix_in_fits(const SpMat &mat, const string &fitsName)
{
    if (mat.rows()*mat.cols() > 2e8)
    {
        LOGLS_WARN(_log, "write_sparse_matrix_in_fits: yout matrix is too large. " << fitsName << " not generated");
        return;
    }
    Mat m(mat.rows(), mat.cols());
    for (int k = 0; k < mat.outerSize(); ++k)
        for (SpMat::InnerIterator it(mat, k); it; ++it)
        {
            m (it.row(), it.col()) = it.value();
        }
    m.writeFits(fitsName);
}

static void write_vect_in_fits(const Eigen::VectorXd &vectorXd, const string &fitsName)
{
    Vect v(vectorXd.size());
    for (int k = 0; k < vectorXd.size(); ++k) v(k) = V(k);
    Mat(v).writeFits(fitsName);
}

#endif


unsigned AstromFit::minimize(const std::string &whatToFit, const double nSigRejCut)
{
    assignIndices(whatToFit);

    // return code can take 3 values :
    // 0 : fit has converged - no more outliers
    // 1 : still some ouliers but chi2 increases
    // 2 : factorization failed
    unsigned returnCode = 0;

    // TODO : write a guesser for the number of triplets
    unsigned nTrip = (_LastNTrip) ? _LastNTrip : 1e6;
    TripletList tList(nTrip);
    Eigen::VectorXd grad(_nParTot);  grad.setZero();

    //Fill the triplets
    LSDerivatives(tList, grad);
    _LastNTrip = tList.size();

    LOGLS_DEBUG(_log, "End of triplet filling, ntrip = " << tList.size());

    SpMat hessian;
    {
#if (TRIPLET_INTERNAL_COORD == COL)
        SpMat jacobian(_nParTot, tList.NextFreeIndex());
        jacobian.setFromTriplets(tList.begin(), tList.end());
        // release memory shrink_to_fit is C++11
        tList.clear(); //tList.shrink_to_fit();
        hessian = jacobian*jacobian.transpose();
#else
        SpMat jacobian(tList.NextRank(), _nParTot);
        jacobian.setFromTriplets(tList.begin(), tList.end());
        // release memory shrink_to_fit is C++11
        tList.clear(); //tList.shrink_to_fit();
        LOGLS_DEBUG(_log, " starting H=JtJ ");
        hessian = jacobian.transpose()*jacobian;
#endif
    }// release the Jacobian


    LOGLS_DEBUG(_log, "Starting factorization, hessian: dim=" << hessian.rows()
               << " nnz=" << hessian.nonZeros()
               << " filling-frac = " << hessian.nonZeros()/sqr(hessian.rows()));

    CholmodSimplicialLDLT2<SpMat> chol(hessian);
    if (chol.info() != Eigen::Success)
    {
        LOGLS_ERROR(_log, "minimize: factorization failed ");
        return 2;
    }

    unsigned tot_outliers = 0;
    double old_chi2 = computeChi2().chi2;

    while (true)
    {
        Eigen::VectorXd delta = chol.solve(grad);
        offsetParams(delta);
        Chi2 current_chi2(computeChi2());
        LOGLS_DEBUG(_log, current_chi2);
        if (current_chi2.chi2 > old_chi2)
        {
            LOGL_WARN(_log, "chi2 went up, exiting outlier rejection loop");
            returnCode = 1;
            break;
        }
        old_chi2 = current_chi2.chi2;

        if (nSigRejCut == 0) break;
        MeasuredStarList moutliers;
        FittedStarList foutliers;
        int n_outliers = findOutliers(nSigRejCut, moutliers, foutliers);
        tot_outliers += n_outliers;
        if (n_outliers == 0) break;
        TripletList tList(1000); // initial allocation size.
        grad.setZero(); // recycle the gradient
        // compute the contributions of outliers to derivatives
        outliersContributions(moutliers, foutliers, tList, grad);
        // actually discard them
        removeMeasOutliers(moutliers);
        removeRefOutliers(foutliers);
        // convert triplet list to eigen internal format
        SpMat h(_nParTot, tList.NextFreeIndex());
        h.setFromTriplets(tList.begin(), tList.end());
        int update_status = chol.update(h, false /* means downdate */);
        LOGLS_DEBUG(_log, "cholmod update_status " << update_status);
        /* The contribution of outliers to the gradient is the opposite
        of the contribution of all other terms, because they add up
         to 0 */
        grad *= -1;
    }

    LOGLS_INFO(_log, "Total number of outliers " << tot_outliers);

    return returnCode;
}

void AstromFit::checkStuff()
{
#if (0)
    const char *what2fit[] = {"Positions", "Distortions", "Refrac",
                              "Positions Distortions", "Positions Refrac",
                              "Distortions Refrac",
                              "Positions Distortions Refrac"
                             };
#endif
    const char *what2fit[] = {"Positions", "Distortions",
                              "Positions Distortions"
                             };
    // DEBUG
    for (int k = 0; k < sizeof(what2fit)/sizeof(what2fit[0]); ++k)
    {
        assignIndices(what2fit[k]);
        TripletList tList(10000);
        Eigen::VectorXd rhs(_nParTot);  rhs.setZero();
        LSDerivatives(tList, rhs);
        SpMat jacobian(_nParTot, tList.NextFreeIndex());
        jacobian.setFromTriplets(tList.begin(), tList.end());
        SpMat hessian = jacobian*jacobian.transpose();
#ifdef STORAGE
        char name[24];
        sprintf(name, "h%d.fits", k);
        write_sparse_matrix_in_fits(hessian, name);
        sprintf(name, "g%d.fits", k);
        write_vect_in_fits(rhs, name);
#endif
        LOGLS_DEBUG(_log, "npar : " << _nParTot << ' ' << _nParDistortions);
    }
}

void AstromFit::makeResTuple(const std::string &tupleName) const
{
    /* cook-up 2 different file names by inserting something just before
       the dot (if any), and within the actual file name. */
    size_t dot = tupleName.rfind('.');
    size_t slash = tupleName.rfind('/');
    if (dot == string::npos || (slash != string::npos && dot < slash))
        dot = tupleName.size();
    std::string meas_tuple(tupleName);
    meas_tuple.insert(dot, "-meas");
    makeMeasResTuple(meas_tuple);
    std::string ref_tuple(tupleName);
    ref_tuple.insert(dot, "-ref");
    makeRefResTuple(ref_tuple);
}

void AstromFit::makeMeasResTuple(const std::string &tupleName) const
{
    std::ofstream tuple(tupleName.c_str());
    tuple << "#xccd: coordinate in CCD" << endl
          << "#yccd: " << endl
          << "#rx:   residual in degrees in TP" << endl
          << "#ry:" << endl
          << "#xtp: transformed coordinate in TP " << endl
          << "#ytp:" << endl
          << "#mag: rough mag" << endl
          << "#jd: Julian date of the measurement" << endl
          << "#rvx: transformed measurement uncertainty " << endl
          << "#rvy:" << endl
          << "#rvxy:" << endl
          << "#color : " << endl
          << "#fsindex: some unique index of the object" << endl
          << "#ra: pos of fitted star" << endl
          << "#dec: pos of fitted star" << endl
          << "#chi2: contribution to Chi2 (2D dofs)" << endl
          << "#nm: number of measurements of this FittedStar" << endl
          << "#chip: chip number" << endl
          << "#visit: visit id" << endl
          << "#end" << endl;
    const CcdImageList &L = _assoc.getCcdImageList();
    for (auto const &i: L)
    {
        const CcdImage &im = *i;
        const MeasuredStarList &cat = im.getCatalogForFit();
        const Mapping *mapping = _distortionModel->getMapping(im);
        const Point &refractionVector = im.RefractionVector();
        double mjd = im.getMjd() - _JDRef;
        for (auto const &is: cat)
        {
            const MeasuredStar &ms = *is;
            if (!ms.IsValid()) continue;
            FatPoint tpPos;
            FatPoint inPos = ms;
            tweakAstromMeasurementErrors(inPos, ms, _posError);
            mapping->TransformPosAndErrors(inPos, tpPos);
            const Gtransfo* sky2TP = _distortionModel->sky2TP(im);
            const FittedStar *fs = ms.GetFittedStar();

            Point fittedStarInTP = transformFittedStar(*fs, sky2TP,
                                   refractionVector,
                                   _refractionCoefficient,
                                   mjd);
            Point res = tpPos - fittedStarInTP;
            double det = tpPos.vx*tpPos.vy - sqr(tpPos.vxy);
            double wxx = tpPos.vy/det;
            double wyy = tpPos.vx/det;
            double wxy = -tpPos.vxy/det;
            //      double chi2 = rx*(wxx*rx+wxy*ry)+ry*(wxy*rx+wyy*ry);
            double chi2 = wxx*res.x*res.x + wyy*res.y*res.y + 2*wxy*res.x*res.y;
            tuple << std::setprecision(9);
            tuple << ms.x << ' ' << ms.y << ' '
                  << res.x << ' ' << res.y << ' '
                  << tpPos.x << ' ' << tpPos.y << ' '
                  << fs->Mag() << ' ' << mjd << ' '
                  << tpPos.vx << ' ' << tpPos.vy << ' ' << tpPos.vxy << ' '
                  << fs->color << ' '
                  << fs->IndexInMatrix() << ' '
                  << fs->x << ' ' << fs->y << ' '
                  << chi2 << ' '
                  << fs->MeasurementCount() << ' '
                  << im.getCcdId() << ' ' << im.getVisit() << endl;
        }// loop on measurements in image
    }// loop on images

}

void AstromFit::makeRefResTuple(const std::string &tupleName) const
{
    std::ofstream tuple(tupleName.c_str());
    tuple << "#ra: coordinates of FittedStar" << endl
          << "#dec: " << endl
          << "#rx:   residual in degrees in TP" << endl
          << "#ry:" << endl
          << "#mag: mag" << endl
          << "#rvx: transformed measurement uncertainty " << endl
          << "#rvy:" << endl
          << "#rvxy:" << endl
          << "#color : " << endl
          << "#fsindex: some unique index of the object" << endl
          << "#chi2: contribution to Chi2 (2D dofs)" << endl
          << "#nm: number of measurements of this FittedStar" << endl
          << "#end" << endl;
    // The following loop is heavily inspired from AstromFit::computeChi2()
    const FittedStarList &fsl = _assoc.fittedStarList;
    TanRaDec2Pix proj(GtransfoLin(), Point(0., 0.));
    for (auto const &i: fsl)
    {
        const FittedStar &fs = *i;
        const RefStar *rs = fs.getRefStar();
        if (rs == nullptr) continue;
        proj.SetTangentPoint(fs);
        // fs projects to (0,0), no need to compute its transform.
        FatPoint rsProj;
        proj.TransformPosAndErrors(*rs, rsProj);
        double rx = rsProj.x; // -fsProj.x (which is 0)
        double ry = rsProj.y;
        double det = rsProj.vx*rsProj.vy - sqr(rsProj.vxy);
        double wxx = rsProj.vy/det;
        double wyy = rsProj.vx/det;
        double wxy = -rsProj.vxy/det;
        double chi2 = wxx*sqr(rx) + 2*wxy*rx*ry + wyy*sqr(ry);
        tuple << std::setprecision(9);
        tuple << fs.x << ' ' << fs.y << ' '
              << rx << ' ' << ry << ' '
              << fs.Mag() << ' '
              << rsProj.vx << ' ' << rsProj.vy << ' ' << rsProj.vxy << ' '
              << fs.color << ' '
              << fs.IndexInMatrix() << ' '
              << chi2 << ' '
              << fs.MeasurementCount() << endl;
    }// loop on FittedStars
}


}
} // end of namespaces
