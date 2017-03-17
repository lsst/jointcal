#ifndef PHOTOMMODEL__H
#define PHOTOMMODEL__H

#include "lsst/jointcal/Eigenstuff.h"
#include <string>
#include <vector>

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;
class MeasuredStar;

//! Interface class for PhotometryFit
class PhotometryModel
{
public :

  /**
   * Assign indices to parameters involved in mappings, starting at firstIndex.
   *
   * @param[in]  whatToFit   String containing parameters to fit.
   * @param[in]  firstIndex  Index to start assigning at.
   *
   * @return     The highest assigned index.
   */
  virtual unsigned assignIndices(const std::string &whatToFit, unsigned firstIndex) = 0;

  /**
   * Offset the parameters by the provided amounts.
   *
   * The shifts are applied according to the indices given in assignIndices.
   *
   * @param[in]  delta  vector of offsets to apply
   */
  virtual void offsetParams(const Eigen::VectorXd &delta) = 0;

  /**
   * Return the "photometric factor" at a given location on a ccdImage.
   *
   * Multiply this by a Calib's flux/magnitude zero-point to get the updated fluxMag0.
   *
   * @param[in]  ccdImage  The ccdImage to get the photometric factor for.
   * @param[in]  where     Possition on ccdImage in ccd coordinates.
   *
   * @return     The photometric factor at the given location on ccdImage.
   */
  virtual double photomFactor(const CcdImage& ccdImage, const Point &where) const =0;

  //! number of parameters to be read in indices.size()
  virtual void getIndicesAndDerivatives(const MeasuredStar &measuredStar,
                                        const CcdImage &ccdImage,
                                        std::vector<unsigned> &indices,
                                        Eigen::VectorXd &D) = 0;

  virtual ~PhotometryModel() {};

};


}} // end of namespaces

#endif /*DISTORTIONSMODEL__H */
