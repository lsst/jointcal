// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H

#include <map>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;

//! Photometric response model which has a single photometric factor per CcdImage.
/*! It considers a full exposure as reference. */
class SimplePhotometryModel : public PhotometryModel {
public:
    SimplePhotometryModel(CcdImageList const &ccdImageList);

    /// No copy or move: there is only ever one instance of a given model (i.e. per ccd+visit)
    SimplePhotometryModel(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel(SimplePhotometryModel &&) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel &&) = delete;

    /**
     * Assign indices to parameters involved in mappings, starting at firstIndex.
     *
     * @param[in]  whatToFit   Ignored.
     * @param[in]  firstIndex  Index to start assigning at.
     *
     * @return     The highest assigned index.
     */
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) override;

    /**
     * Offset the parameters by the provided amounts.
     *
     * The shifts are applied according to the indices given in AssignIndices.a
     *
     * @param[in]  delta  vector of offsets to apply
     */
    void offsetParams(Eigen::VectorXd const &delta) override;

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
    double photomFactor(CcdImage const &ccdImage, Point const &where = Point()) const override;

    void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) override;

    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) override;

private:
    typedef std::map<CcdImage const *, std::unique_ptr<PhotometryMapping>> MapType;
    MapType _myMap;

    /// Return the mapping associated with this ccdImage. name is a descriptor for error messages.
    PhotometryMapping *findMapping(CcdImage const &ccdImage, std::string name) const override;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
