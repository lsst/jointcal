#ifndef SIMPLEPHOTOMMODEL__H
#define SIMPLEPHOTOMMODEL__H

#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/Point.h"
#include <map>

namespace lsst {
namespace jointcal {

class CcdImageList;
class CcdImage;
class Point;

//! Photometric response model which has a single photometric factor per CcdImage.
/*! It considers a full exposure as reference. */
 class SimplePhotometryModel : public PhotometryModel
{

  struct PhotomStuff
  {
    unsigned index;
    double factor;
    bool fixed;
    PhotomStuff(const unsigned i=0, const double f=1) : index(i), factor(f), fixed(false) {};
  };

  typedef std::map<const CcdImage*,PhotomStuff> mapType;
  mapType _myMap;

  PhotomStuff& find(const CcdImage &ccdImage);
  const PhotomStuff& find(const CcdImage &ccdImage) const;

public :

  SimplePhotometryModel(const CcdImageList &ccdImageList);

  /**
   * Assign indices to parameters involved in mappings, starting at firstIndex.
   *
   * @param[in]  whatToFit   Ignored.
   * @param[in]  firstIndex  Index to start assigning at.
   *
   * @return     The highest assigned index.
   */
  unsigned assignIndices(const std::string &whatToFit, unsigned firstIndex);

  /**
   * Offset the parameters by the provided amounts.
   *
   * The shifts are applied according to the indices given in AssignIndices.a
   *
   * @param[in]  delta  vector of offsets to apply
   */
  void offsetParams(const Eigen::VectorXd &delta);

  /**
   * Return the "photometric factor" for this ccdImage.
   *
   * Multiply this by a Calib's flux/magnitude zero-point to get the updated fluxMag0.
   *
   * @param[in]  ccdImage  The ccdImage to get the photometric factor for.
   * @param[in]  where     Ignored
   *
   * @return     The photometric factor at the given location on ccdImage.
   */
  double photomFactor(const CcdImage& ccdImage, const Point &where=Point()) const;

  void getIndicesAndDerivatives(const MeasuredStar &measuredStar,
                                const CcdImage &ccdImage,
                                std::vector<unsigned> &indices,
                                Eigen::VectorXd &D);

};


}} // end of namespaces

#endif /*SIMPLEPHOTOMMODEL__H */
