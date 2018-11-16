// -*- LSST-C++ -*-
/*
 * This file is part of jointcal.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "lsst/log/Log.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/PhotometryModel.h"

namespace lsst {
namespace jointcal {

bool PhotometryModel::validate(CcdImageList const &ccdImageList) const {
    bool check = true;
    for (auto const &ccdImage : ccdImageList) {
        // Don't short circuit so that we log every place the model is negative.
        check &= checkPositiveOnBBox(*ccdImage);
    }
    return check;
}

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
