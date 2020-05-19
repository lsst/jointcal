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
        # We don't want the absolute astrometry to become significantly worse
        # than the single-epoch astrometry (about 0.040").
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 49e-3*u.arcsecond

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
            dist_rms_relative = 11e-3*u.arcsecond
            pa1 = 0.014
            metrics = {'collected_astrometry_refStars': 1770,
                       'collected_photometry_refStars': 1770,
                       'selected_astrometry_refStars': 747,
                       'selected_photometry_refStars': 747,
                       'associated_astrometry_fittedStars': 2269,
                       'associated_photometry_fittedStars': 2269,
                       'selected_astrometry_fittedStars': 1408,
                       'selected_photometry_fittedStars': 1408,
                       'selected_astrometry_ccdImages': 12,
                       'selected_photometry_ccdImages': 12,
                       'astrometry_final_chi2': 1609.29,
                       'astrometry_final_ndof': 3332,
                       'photometry_final_chi2': 3632.26,
                       'photometry_final_ndof': 1693
                       }

            self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

            # Check for the existence of the chi2 contribution files.
            expected = ['photometry_initial_chi2-0_r', 'astrometry_initial_chi2-0_r',
                        'photometry_final_chi2-0_r', 'astrometry_final_chi2-0_r']
            for partial in expected:
                name = os.path.join(tempdir, partial+'-ref.csv')
                self.assertTrue(os.path.exists(name), msg="Did not find file %s"%name)
                name = os.path.join(tempdir, partial+'-meas.csv')
                self.assertTrue(os.path.exists(name), msg='Did not find file %s'%name)

            expected = ["initialAstrometryModel.txt", "initialPhotometryModel.txt"]
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
        dist_rms_relative = 12e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 1770,
                   'selected_astrometry_refStars': 747,
                   'associated_astrometry_fittedStars': 2269,
                   'selected_astrometry_fittedStars': 1408,
                   'selected_astrometry_ccdImages': 12,
                   'astrometry_final_chi2': 1714.8,
                   'astrometry_final_ndof': 3434,
                   }

        return dist_rms_relative, metrics

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.writeInitialModel = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            self.config.debugOutputPath = tempdir

            self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)
            filename = os.path.join(tempdir, "initialAstrometryModel.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedAstrometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect astrometry.
        """
        relative_error, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.astrometryDoRankUpdate = False

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_4sigma_outliers(self):
        """4 sigma outlier rejection means fewer available sources after the
        fitter converges, resulting in a smaller ndof and chi2.
        """
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self.config.outlierRejectSigma = 4
        metrics['astrometry_final_chi2'] = 1288.64
        metrics['astrometry_final_ndof'] = 3232

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryReferenceUncertainty_smaller(self):
        """Test with a smaller fake reference uncertainty: chi2 will be higher."""
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/astrometryReferenceErr-config.py')
        self.configfiles.append(test_config)
        metrics['astrometry_final_chi2'] = 11522.9
        metrics['astrometry_final_ndof'] = 3406

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
        pa1 = 0.017
        metrics = {'collected_photometry_refStars': 1770,
                   'selected_photometry_refStars': 747,
                   'associated_photometry_fittedStars': 2269,
                   'selected_photometry_fittedStars': 1408,
                   'selected_photometry_ccdImages': 12,
                   'photometry_final_chi2': 3292.08,
                   'photometry_final_ndof': 1622
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.writeInitialModel = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            self.config.debugOutputPath = tempdir

            self._testJointcalTask(2, None, None, pa1, metrics=metrics)
            filename = os.path.join(tempdir, "initialPhotometryModel.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedPhotometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect photometry.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryDoRankUpdate = False

        # The constrainedPhotometry model is not purely linear, so a small
        # change in final chi2 is possible.
        metrics['photometry_final_chi2'] = 3297.34

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
        metrics['photometry_final_chi2'] = 3324.2
        metrics['photometry_final_ndof'] = 1625

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
        metrics['selected_photometry_refStars'] = 214
        metrics['associated_photometry_fittedStars'] = 270
        metrics['selected_photometry_fittedStars'] = 245
        metrics['photometry_final_chi2'] = 373.141
        metrics['photometry_final_ndof'] = 254

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedMagnitude_no_astrometry(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryModel = "constrainedMagnitude"

        # The resulting fit should be close to the constrainedFlux model:
        # there are few CCDs and 2 visits, so there's not a lot of complexity
        # in this case to distinguish the flux vs. magnitude models.
        metrics['photometry_final_chi2'] = 3332.76

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedFlux_pedestal(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryErrorPedestal = 0.02

        # We're allowing more error in the fit, so PA1 may be worse.
        pa1 = 0.021
        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 2246.58
        # ndof shouldn't change much; slightly different likelihood contours
        metrics['photometry_final_ndof'] = 1624

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedMagnitude_pedestal(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryModel = "constrainedMagnitude"
        self.config.photometryErrorPedestal = 0.02

        # We're allowing more error in the fit, so PA1 may be worse.
        pa1 = 0.024
        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 2243.56
        # ndof shouldn't change much; slightly different likelihood contours
        metrics['photometry_final_ndof'] = 1617

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
