#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>

#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/AstroUtils.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/jointcal/Gtransfo.h"

namespace lsst {
namespace jointcal {

Frame applyTransfo(const Frame& inputframe, const Gtransfo& gtransfo, const WhichTransformed which) {
    // 2 opposite corners
    double xtmin1, xtmax1, ytmin1, ytmax1;
    gtransfo.apply(inputframe.xMin, inputframe.yMin, xtmin1, ytmin1);
    gtransfo.apply(inputframe.xMax, inputframe.yMax, xtmax1, ytmax1);
    Frame fr1(std::min(xtmin1, xtmax1), std::min(ytmin1, ytmax1), std::max(xtmin1, xtmax1),
              std::max(ytmin1, ytmax1));
    // 2 other corners
    double xtmin2, xtmax2, ytmin2, ytmax2;
    gtransfo.apply(inputframe.xMin, inputframe.yMax, xtmin2, ytmax2);
    gtransfo.apply(inputframe.xMax, inputframe.yMin, xtmax2, ytmin2);
    Frame fr2(std::min(xtmin2, xtmax2), std::min(ytmin2, ytmax2), std::max(xtmin2, xtmax2),
              std::max(ytmin2, ytmax2));

    if (which == SmallFrame) return fr1 * fr2;
    return fr1 + fr2;
}
}  // namespace jointcal
}  // namespace lsst
