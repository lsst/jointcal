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

import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestDECAM(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 63e-3*u.arcsecond

        do_plot = False

        # center of the decam validation_data catalog
        center = lsst.afw.geom.SpherePoint(150.1191666, 2.20583333, lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'decam')
        all_visits = [176837, 176846]
        # Skipping ccd=13,14 because they mangle the results (14 is missing in one visit, 13 is bad).
        ccdnums = '^'.join(str(x) for x in (10, 11, 12, 15, 16, 17, 18))
        other_args = ['ccdnum=' + ccdnums, ]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot,
                        log_level="DEBUG")

        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/decam-config.py')
        self.configfiles.append(test_config)

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")

        # See Readme for an explanation of these empirical values.
        # NOTE: the photometry and astrometry refstars numbers are different
        # here because the SDSS catalogs have some sources with bad fluxes;
        # those are rejected for photometric calibration, but not astrometric.
        relative_error = 19e-3*u.arcsecond
        pa1 = 0.14
        metrics = {'collected_astrometry_refStars': 8869,
                   'collected_photometry_refStars': 8858,
                   'selected_astrometry_refStars': 1604,
                   'selected_photometry_refStars': 1604,
                   'associated_astrometry_fittedStars': 6749,
                   'associated_photometry_fittedStars': 6749,
                   'selected_astrometry_fittedStars': 2709,
                   'selected_photometry_fittedStars': 2709,
                   'selected_astrometry_ccdImages': 14,
                   'selected_photometry_ccdImages': 14,
                   'astrometry_final_chi2': 2838.0,
                   'astrometry_final_ndof': 5562,
                   'photometry_final_chi2': 4810.5,
                   'photometry_final_ndof': 2910,
                   }

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        """Set default values for the constrainedAstrometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "constrained"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        relative_error = 17e-3*u.arcsecond
        pa1 = None
        metrics = {'collected_astrometry_refStars': 8869,
                   'selected_astrometry_refStars': 1604,
                   'associated_astrometry_fittedStars': 6749,
                   'selected_astrometry_fittedStars': 2709,
                   'selected_astrometry_ccdImages': 14,
                   'astrometry_final_chi2': 3040.84,
                   'astrometry_final_ndof': 5748,
                   }
        return relative_error, pa1, metrics

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        relative_error, pa1, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry_no_photometry()

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't significantly affect astrometry.
        """
        relative_error, pa1, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry_no_photometry()
        self.config.astrometryDoRankUpdate = False

        metrics['astrometry_final_chi2'] = 3015.61
        metrics['astrometry_final_ndof'] = 5738

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    @unittest.skip("DM-14439 : This test produces different chi2/ndof on centos6 and centos7/macOS.")
    def test_jointcalTask_2_visits_constrainedAstrometry_4sigma_outliers(self):
        """4 sigma outlier rejection means a fewer available sources, resulting
        in a smaller ndof and chi2.
        """
        relative_error, pa1, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry_no_photometry()
        self.config.outlierRejectSigma = 4
        metrics['astrometry_final_chi2'] = 1506.7
        metrics['astrometry_final_ndof'] = 3682

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        """Set default values for the constrainedPhotometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "constrainedFlux"
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        self.other_args.extend(['--config', 'writeChi2ContributionFiles=True'])

        # See Readme for an explanation of these empirical values.
        pa1 = 0.11
        metrics = {'collected_photometry_refStars': 4865,
                   'selected_photometry_refStars': 661,
                   'associated_photometry_fittedStars': 6749,
                   'selected_photometry_fittedStars': 2044,
                   'selected_photometry_ccdImages': 14,
                   'photometry_final_chi2': 3066.22,
                   'photometry_final_ndof': 1998,
                   }

        return pa1, metrics

    @unittest.skip("DM-14439 : This test produces different chi2/ndof on Linux and macOS.")
    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry_fluxErr_0(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry()

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    @unittest.skip("DM-14439 : This test produces different chi2/ndof on Linux and macOS.")
    def test_jointcalTask_2_visits_constrainedPhotometry_flagged_selector(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry()
        self.config.sourceSelector.name = 'flagged'
        # Reduce warnings due to flaggedSourceSelector having fewer sources than astrometrySourceSelector.
        self.config.minMeasuredStarsPerCcd = 40

        # See Readme for an explanation of these empirical values.
        metrics = {'collected_photometry_refStars': 4865,
                   'selected_photometry_refStars': 551,
                   'associated_photometry_fittedStars': 860,
                   'selected_photometry_fittedStars': 593,
                   'selected_photometry_ccdImages': 14,
                   'photometry_final_chi2': 817.124,
                   'photometry_final_ndof': 607,
                   }

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    @unittest.skip("DM-14439 : This test produces different chi2/ndof on Linux and macOS.")
    def test_jointcalTask_2_visits_constrainedPhotometry_no_rank_update(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry()
        self.config.photometryDoRankUpdate = False

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    @unittest.skip("DM-14439 : This test produces different chi2/ndof on Linux and macOS.")
    def test_jointcalTask_2_visits_constrainedPhotometry_fluxErr_3(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry()
        self.config.fluxError = 0.03
        metrics['photometry_final_chi2'] = 1818.92
        metrics['photometry_final_ndof'] = 2010

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
