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

import unittest
import os

from astropy import units as u

import lsst.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestDECAM(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    """NOTE: the data used in this test has some bad input astrometry
    (particularly in visit 176837) so the data coming in to jointcal is of
    poor quality. The fitter logs many warnings about not having enough
    sources, and not being able to produce a good inverse model for the
    WCS output. This data does not represent a good test case for
    jointcal's fitter, but may serve as a useful test case for how to fail
    on a bad fit.
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_decam')
        except LookupError:
            raise unittest.SkipTest("obs_decam not setup")

    def setUp(self):
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 63e-3*u.arcsecond

        do_plot = False

        # center of the decam validation_data catalog
        center = lsst.geom.SpherePoint(150.1191666, 2.20583333, lsst.geom.degrees)
        radius = 3*lsst.geom.degrees

        input_dir = os.path.join(self.data_dir, 'decam')
        all_visits = [176837, 176846]
        ccdnums = '^'.join(str(x) for x in (10, 11, 12, 14, 15, 16, 17, 18))
        other_args = ['ccdnum=' + ccdnums, ]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot,
                        log_level="DEBUG")

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"

        # See Readme for an explanation of these empirical values.
        relative_error = 19e-3*u.arcsecond
        pa1 = 0.24
        metrics = {'collected_astrometry_refStars': 7403,
                   'collected_photometry_refStars': 149481,
                   'selected_astrometry_refStars': 467,
                   'selected_photometry_refStars': 4632,
                   'associated_astrometry_fittedStars': 7839,
                   'associated_photometry_fittedStars': 7839,
                   'selected_astrometry_fittedStars': 702,
                   'selected_photometry_fittedStars': 4637,
                   'selected_astrometry_ccdImages': 15,
                   'selected_photometry_ccdImages': 15,
                   'astrometry_final_chi2': 1.21494,
                   'astrometry_final_ndof': 560,
                   'photometry_final_chi2': 16529,
                   'photometry_final_ndof': 3634,
                   }

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
