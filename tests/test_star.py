# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import warnings

import numpy as np
import astropy.units as u
import astropy.time
from astropy.coordinates import SkyCoord

import unittest
import lsst.utils.tests

import lsst.jointcal.star


class TestProperMotion(lsst.utils.tests.TestCase):
    """Tests of applying proper motion correction to stars."""
    def setUp(self):
        # These tests use times that are outside what ERFA considers "reasonable";
        # ignore any warnings it emits about times and distances.
        warnings.filterwarnings("ignore", module="erfa")
        self.ra = 300 * u.degree
        self.dec = 80 * u.degree
        self.pm_ra = 1000 * u.mas / u.yr  # ra*cos(dec)
        self.pm_dec = 2000 * u.mas / u.yr
        flux = 10  # need something for the constructor, but value is irrelevant here
        # 100 year baseline
        self.baseline_epoch = astropy.time.Time("2000-01-01 00:00:05", scale="tai")
        self.observed_epoch = astropy.time.Time("2100-01-01 00:00:05", scale="tai")
        self.dt = (self.observed_epoch.tai - self.baseline_epoch.tai).to(astropy.units.yr)
        self.coord = SkyCoord(self.ra, self.dec, frame="icrs",
                              pm_ra_cosdec=self.pm_ra, pm_dec=self.pm_dec,
                              obstime=self.observed_epoch)
        self.star = lsst.jointcal.star.BaseStar(self.ra.value, self.dec.value,
                                                flux, flux*0.001)
        self.refStar = lsst.jointcal.star.RefStar(self.ra.value, self.dec.value,
                                                  flux, flux*0.001)
        self.star.vx = self.ra.to_value(u.degree) * 0.01
        self.star.vy = self.dec.to_value(u.degree) * 0.01
        self.properMotion = lsst.jointcal.star.ProperMotion(self.pm_ra.to_value(u.radian/u.yr),
                                                            self.pm_dec.to_value(u.radian/u.yr),
                                                            self.pm_ra.to_value(u.radian/u.yr)*0.01,
                                                            self.pm_dec.to_value(u.radian/u.yr)*0.01)

        # Test points on the whole sphere, all with the same large proper motion value.
        np.random.seed(100)
        n = 100
        ras = np.random.random(n)*360 * u.degree
        # uniformly distributed in declination
        decs = np.arccos(np.random.random(n)*2 - 1) * u.degree
        pm_ra = np.ones(n) * 300 * u.mas / u.yr
        pm_dec = np.ones(n) * 400 * u.mas / u.yr
        self.coords = SkyCoord(ras, decs, frame="icrs", pm_ra_cosdec=pm_ra, pm_dec=pm_dec,
                               obstime=self.observed_epoch)

    def test_apply_no_proper_motion(self):
        """If refCat as no ProperMotion set, applyProperMotion() should change
        nothing.
        """
        result = self.refStar.applyProperMotion(self.star, self.dt.value)
        self.assertEqual(result.x, self.star.x)
        self.assertEqual(result.y, self.star.y)
        # TODO? astropy SkyCoord does not include coordinate errors, or error propogation.
        # How do I test it?
        # self.assertEqual(result.vx, self.star.vx)
        # self.assertEqual(result.vy, self.star.vy)

    def test_apply_one(self):
        """Test apply on a single coordinate (useful for debugging).
        """
        expect = self.coord.apply_space_motion(dt=self.dt)

        self.refStar.setProperMotion(self.properMotion)
        result = self.refStar.applyProperMotion(self.star, self.dt.value)
        # original star should not be changed:
        self.assertEqual(self.ra.to_value(u.degree), self.star.x)
        self.assertEqual(self.dec.to_value(u.degree), self.star.y)
        self.assertEqual(self.ra.to_value(u.degree)*0.01, self.star.vx)
        self.assertEqual(self.dec.to_value(u.degree)*0.01, self.star.vy)

        # 1e-9 deg == 3.6 microarcsec; that's pretty good accuracy over a 100 year baseline.
        self.assertFloatsAlmostEqual(result.x, expect.ra.to_value(u.degree), rtol=1e-9)
        self.assertFloatsAlmostEqual(result.y, expect.dec.to_value(u.degree), rtol=1e-9)
        # TODO? astropy SkyCoord does not include coordinate errors, or error propogation.
        # How do I test it?
        # self.assertEqual(result.vx, expect.vx)
        # self.assertEqual(result.vy, expect.vy)

    def test_apply_many(self):
        """Test apply over a range of points on the sphere.
        """
        expect = self.coords.apply_space_motion(dt=self.dt)

        ras = np.zeros(len(expect))
        decs = np.zeros(len(expect))
        for i, x in enumerate(self.coords):
            self.coords
            star = lsst.jointcal.star.BaseStar(x.ra.value, x.dec.value,
                                               100, 100*0.001)
            refStar = lsst.jointcal.star.RefStar(x.ra.value, x.dec.value,
                                                 100, 100*0.001)
            star.vx = x.ra.to_value(u.degree) * 0.01
            star.vy = x.dec.to_value(u.degree) * 0.01
            properMotion = lsst.jointcal.star.ProperMotion(x.pm_ra_cosdec.to_value(u.radian/u.yr),
                                                           x.pm_dec.to_value(u.radian/u.yr),
                                                           x.pm_ra_cosdec.to_value(u.radian/u.yr)*0.01,
                                                           x.pm_dec.to_value(u.radian/u.yr)*0.01)
            refStar.setProperMotion(properMotion)
            result = refStar.applyProperMotion(star, self.dt.value)
            ras[i] = result.x
            decs[i] = result.y
        self.assertFloatsAlmostEqual(ras, expect.ra.to_value(u.degree), rtol=1e-7)
        self.assertFloatsAlmostEqual(decs, expect.dec.to_value(u.degree), rtol=6e-7)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
