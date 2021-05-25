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
import tempfile

from astropy import units as u

import lsst.geom
import lsst.utils
import lsst.pex.exceptions
import lsst.pex.config

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHT(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_cfht')
        except LookupError:
            raise unittest.SkipTest("obs_cfht not setup")

    def setUp(self):
        # NOTE: refcat-comparison RMS error is worse now, because the
        # comparison code is not applying the proper motion data.
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 70e-3*u.arcsecond

        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.geom.SpherePoint(214.884832, 52.6622199, lsst.geom.degrees)
        radius = 3*lsst.geom.degrees

        input_dir = os.path.join(self.data_dir, 'cfht')
        all_visits = [849375, 850587]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot,
                        log_level="DEBUG")

    def test_jointcalTask_2_visits(self):
        """Test the simple models with two visits and check that some debug
        output files also get created.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        self.config.writeInitialModel = True  # write the initial models
        # to test whether we got the expected chi2 contribution files.
        self.config.writeChi2FilesInitialFinal = True
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            self.config.debugOutputPath = tempdir

            # See Readme for an explanation of these empirical values.
            dist_rms_relative = 16e-3*u.arcsecond
            pa1 = 0.049
            metrics = {'collected_astrometry_refStars': 867,
                       'collected_photometry_refStars': 11570,
                       'selected_astrometry_refStars': 332,
                       'selected_photometry_refStars': 2225,
                       'associated_astrometry_fittedStars': 2272,
                       'associated_photometry_fittedStars': 2272,
                       'selected_astrometry_fittedStars': 1229,
                       'selected_photometry_fittedStars': 2232,
                       'selected_astrometry_ccdImages': 12,
                       'selected_photometry_ccdImages': 12,
                       'astrometry_final_chi2': 1159.736,
                       'astrometry_final_ndof': 1872,
                       'photometry_final_chi2': 11624.3,
                       'photometry_final_ndof': 2778
                       }

            self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

            # Check for the existence of the chi2 contribution files.
            expected = ['photometry_initial_chi2-0_r.MP9601', 'astrometry_initial_chi2-0_r.MP9601',
                        'photometry_final_chi2-0_r.MP9601', 'astrometry_final_chi2-0_r.MP9601']
            for partial in expected:
                name = os.path.join(tempdir, partial+'-ref.csv')
                self.assertTrue(os.path.exists(name), msg="Did not find file %s"%name)
                name = os.path.join(tempdir, partial+'-meas.csv')
                self.assertTrue(os.path.exists(name), msg='Did not find file %s'%name)

            expected = ["initial_astrometry_model-0_r.MP9601.txt", "initial_photometry_model-0_r.MP9601.txt"]
            for name in expected:
                fullpath = os.path.join(tempdir, name)
                self.assertTrue(os.path.exists(fullpath), msg=f"Did not find file {fullpath}")

    def setup_jointcalTask_2_visits_constrainedAstrometry(self):
        """Set default values for the constrainedAstrometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "constrained"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 16e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 867,
                   'selected_astrometry_refStars': 332,
                   'associated_astrometry_fittedStars': 2272,
                   'selected_astrometry_fittedStars': 1229,
                   'selected_astrometry_ccdImages': 12,
                   'astrometry_final_chi2': 1124.0575,
                   'astrometry_final_ndof': 1912,
                   }

        return dist_rms_relative, metrics

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.writeInitialModel = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            self.config.debugOutputPath = tempdir

            self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)
            filename = os.path.join(tempdir, "initial_astrometry_model-0_r.MP9601.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedAstrometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect astrometry.
        """
        relative_error, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        metrics['astrometry_final_chi2'] = 1069.348
        metrics['astrometry_final_ndof'] = 1900

        self.config.astrometryDoRankUpdate = False

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_4sigma_outliers(self):
        """4 sigma outlier rejection means fewer available sources after the
        fitter converges, resulting in a smaller ndof and chi2.
        """
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.outlierRejectSigma = 4
        metrics['astrometry_final_chi2'] = 757.027
        metrics['astrometry_final_ndof'] = 1732

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryOutlierRelativeTolerance(self):
        """Test that astrometryOutlierRelativeTolerance changes the fit. Setting
        1% for the astrometryOutlierRelativeTolerance will result in higher chi2
        and ndof.
        """
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.astrometryOutlierRelativeTolerance = 0.01
        metrics['astrometry_final_chi2'] = 1427.11
        metrics['astrometry_final_ndof'] = 1990

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryReferenceUncertainty_smaller(self):
        """Test with a smaller fake reference uncertainty: chi2 will be higher."""
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/astrometryReferenceErr-config.py')
        self.configfiles.append(test_config)
        metrics['astrometry_final_chi2'] = 1275.70
        metrics['astrometry_final_ndof'] = 2062

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryReferenceUncertainty_None_fails(self):
        """The default `None` should fail for the existing refcats that have no position errors."""
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        # This is the default, but we override it in tests/config/config.py,
        # because none of the existing test refcats have errors. So we have to
        # re-override it here.
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/astrometryReferenceErr-None-config.py')
        self.configfiles.append(test_config)
        with self.assertRaises(lsst.pex.config.FieldValidationError):
            self._testJointcalTask(2, None, None, None)

    def setup_jointcalTask_2_visits_constrainedPhotometry(self):
        """Set default values for the constrainedPhotometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "constrainedFlux"
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        # See Readme for an explanation of these empirical values.
        pa1 = 0.09
        metrics = {'collected_photometry_refStars': 11570,
                   'selected_photometry_refStars': 2225,
                   'associated_photometry_fittedStars': 2272,
                   'selected_photometry_fittedStars': 2232,
                   'selected_photometry_ccdImages': 12,
                   'photometry_final_chi2': 11649.7,
                   'photometry_final_ndof': 2729
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.writeInitialModel = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            self.config.debugOutputPath = tempdir

            self._testJointcalTask(2, None, None, pa1, metrics=metrics)
            filename = os.path.join(tempdir, "initial_photometry_model-0_r.MP9601.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedPhotometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect photometry.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryDoRankUpdate = False

        # The constrainedPhotometry model is not purely linear, so a small
        # change in final chi2 is possible.
        metrics['photometry_final_chi2'] = 11396.27
        metrics['photometry_final_ndof'] = 2716

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_lineSearch(self):
        """Activating the line search should only slightly change the chi2.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.allowLineSearch = True

        # Activating line search for constrainedPhotometry should result in
        # nearly the same final fit (the system is somewhat non-linear, so it
        # may not be exactly the same: check the "Line search scale factor"
        # lines in the DEBUG log for values that are not ~1 for proof).
        pa1 = 0.14
        metrics['photometry_final_chi2'] = 10500.4
        metrics['photometry_final_ndof'] = 2714

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_flagged(self):
        """Test the use of the FlaggedSourceSelector."""
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/cfht-flagged-config.py')
        self.configfiles.append(test_config)
        # Reduce warnings due to flaggedSourceSelector having fewer sources than astrometrySourceSelector.
        self.config.minMeasuredStarsPerCcd = 30
        self.config.minRefStarsPerCcd = 20

        pa1 = 0.026
        metrics['selected_photometry_refStars'] = 265
        metrics['associated_photometry_fittedStars'] = 265
        metrics['selected_photometry_fittedStars'] = 265
        metrics['photometry_final_chi2'] = 392.816
        metrics['photometry_final_ndof'] = 294

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedMagnitude_no_astrometry(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryModel = "constrainedMagnitude"

        # The resulting fit should be close to the constrainedFlux model:
        # there are few CCDs and 2 visits, so there's not a lot of complexity
        # in this case to distinguish the flux vs. magnitude models.
        metrics['photometry_final_chi2'] = 9455.06
        metrics['photometry_final_ndof'] = 2710

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedFlux_pedestal(self):
        """Test that forcing a systematic flux error results in a lower chi2.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        # median fluxErr/flux in ps1 is 0.11, so we have to make this bigger
        # than that to actually allow more slop in the fit.
        self.config.photometryErrorPedestal = 0.2

        # We're allowing more error in the fit, so PA1 may be worse.
        pa1 = 0.21
        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 3214.78
        # ndof may change; slightly different likelihood contours, and fewer
        # reference sources rejected.
        metrics['photometry_final_ndof'] = 3154

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedMagnitude_pedestal(self):
        """Test that forcing a systematic flux error results in a lower chi2.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryModel = "constrainedMagnitude"
        # median fluxErr/flux in ps1 is 0.11, so we have to make this bigger
        # than that to actually allow more slop in the fit.
        self.config.photometryErrorPedestal = 0.2

        # We're allowing more error in the fit, so PA1 may be worse.
        pa1 = 0.19
        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 3092.39
        # ndof may change; slightly different likelihood contours, and fewer
        # reference sources rejected.
        metrics['photometry_final_ndof'] = 3129

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
