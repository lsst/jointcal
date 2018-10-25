#include "lsst/log/Log.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotometryModel.h"

namespace lsst {
namespace jointcal {

bool PhotometryModel::checkPositiveOnBBox(CcdImage const &ccdImage) const {
    bool check = true;
    auto bbox = ccdImage.getImageFrame();
    for (auto const &x : {bbox.xMin, bbox.getCenter().x, bbox.xMax})
        for (auto const &y : {bbox.yMin, bbox.getCenter().y, bbox.yMax}) {
            // flux, fluxErr, instFlux, instFluxErr all 1
            jointcal::MeasuredStar star(jointcal::BaseStar(x, y, 1, 1));
            star.setInstFluxAndErr(1, 1);
            double result = transform(ccdImage, star);
            // Don't short circuit so that we log every place the model is negative.
            if (result < 0) {
                LOGLS_WARN(_log, "CcdImage " << ccdImage.getName() << " is negative at (" << x << "," << y
                                             << "): " << result);
                check = false;
            }
        }
    return check;
}

}  // namespace jointcal
}  // namespace lsst
