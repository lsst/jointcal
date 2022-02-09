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

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE test_trans

// The boost unit test header
#include "boost/test/unit_test.hpp"

#include "lsst/jointcal/StarMatch.h"
#include "lsst/jointcal/AstrometryTransform.h"
#include "lsst/jointcal/Frame.h"
#include "lsst/geom/Angle.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/fits.h"
#include "lsst/daf/base.h"

#include <cstdlib> /* for getenv */

// NOTE: turn this flag on to raise exceptions on floating point errors.
// NOTE: this only works on GNU/Linux (fenableexcept is not C++ standard).
// #define DUMP_CORE_ON_FPE
#ifdef DUMP_CORE_ON_FPE
#define _GNU_SOURCE 1
#define __USE_GNU
#include <fenv.h>
static void __attribute__((constructor)) trapfpe() {
    // Enable some exceptions.  At startup all exceptions are masked.
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#endif

namespace jointcal = lsst::jointcal;
namespace afwImg = lsst::afw::image;

/* Test jointcal::TanSipPixelToRaDec::apply against afwImg::Wcs::pixelToSky */

BOOST_AUTO_TEST_SUITE(test_transforms)

BOOST_AUTO_TEST_CASE(test_wcs) {
    std::string fileName = "tests/header_only.fits";

    lsst::afw::fits::Fits file(fileName, "r", 0);
    auto propSet = lsst::afw::fits::readMetadata(fileName);
    auto skyWcs = lsst::afw::geom::makeSkyWcs(*propSet);
    jointcal::AstrometryTransformSkyWcs astrometryTransformSkyWcs(skyWcs);

    jointcal::Point where(100., 200.);
    jointcal::Point outPol = astrometryTransformSkyWcs.apply(where);
    std::cout << std::setprecision(12) << "Poloka : " << outPol.x << ' ' << outPol.y << std::endl;

    lsst::geom::Point2D whereSame(100., 200.);
    auto skyPos = skyWcs->pixelToSky(whereSame);
    lsst::geom::Point2D outDeg = skyPos.getPosition(lsst::geom::degrees);
    std::cout << "Stack : " << outDeg[0] << ' ' << outDeg[1] << std::endl;

    BOOST_CHECK_CLOSE(outPol.x, outDeg[0], .000001);
    BOOST_CHECK_CLOSE(outPol.y, outDeg[1], .000001);
}

/* test the AstrometryTransformPolynomial::fit routine */

BOOST_AUTO_TEST_CASE(test_polyfit) {
    std::string fileName = "tests/header_only.fits";

    lsst::afw::fits::Fits file(fileName, "r", 0);
    auto propSet = lsst::afw::fits::readMetadata(fileName);
    auto skyWcs = lsst::afw::geom::makeSkyWcs(*propSet);
    jointcal::AstrometryTransformSkyWcs astrometryTransformSkyWcs(skyWcs);

    jointcal::StarMatchList sml;
    jointcal::BaseStarList bsl1, bsl2;

    for (double x = 10; x < 2000; x += 120)
        for (double y = 20; y < 4000; y += 160) {
            auto s1 = std::make_shared<jointcal::BaseStar>(x, y, 1, 0.01);
            s1->vx = 0.1;
            s1->vy = 0.2;
            s1->vxy = 0.05;
            auto s2 = std::make_shared<jointcal::BaseStar>();
            astrometryTransformSkyWcs.transformPosAndErrors(*s1, *s2);
            bsl1.push_back(s1);
            bsl2.push_back(s2);
            sml.push_back(jointcal::StarMatch(*s1, *s2, s1, s2));
        }
    jointcal::AstrometryTransformPolynomial pol(3);
    double chi2 = pol.fit(sml);
    std::cout << " chi2/ndf " << chi2 << '/' << sml.size() - pol.getNpar() << std::endl;
    // since there is no noise, the chi2 should be very very small:
    BOOST_CHECK(fabs(chi2) < 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
