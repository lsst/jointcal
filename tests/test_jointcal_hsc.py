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
import lsst.pex.config
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestHSC(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 53e-3*u.arcsecond

        do_plot = False

        # center of the hsc validation_data catalog
        center = lsst.afw.geom.SpherePoint(320.367492, 0.3131554, lsst.afw.geom.degrees)
        radius = 5*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'hsc')
        all_visits = [903334, 903336, 903338, 903342, 903344, 903346, 903986, 903988, 903990, 904010, 904014]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-config.py')
        self.configfiles.append(test_config)

    def test_jointcalTask_2_visits_simple(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.024
        metrics = {'collected_astrometry_refStars': 2948,
                   'collected_photometry_refStars': 2948,
                   'selected_astrometry_refStars': 657,
                   'selected_photometry_refStars': 657,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 851,
                   'selected_photometry_fittedStars': 851,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_final_chi2': 819.07,
                   'astrometry_final_ndof': 2134,
                   'photometry_final_chi2': 1811.45,
                   'photometry_final_ndof': 1072
                   }
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_11_visits_simple_astrometry_no_photometry(self):
        """Test 11 visits with different filters.
        Testing photometry doesn't make sense for this currently.
        """

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = None  # pa1 = 0.134
        metrics = {'collected_astrometry_refStars': 4926,
                   'selected_astrometry_refStars': 1346,
                   'associated_astrometry_fittedStars': 2908,
                   'selected_astrometry_fittedStars': 2276,
                   'selected_astrometry_ccdImages': 33,
                   'astrometry_final_chi2': 8207.62,
                   'astrometry_final_ndof': 14856,
                   }
        self._testJointcalTask(11, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_simplePhotometry(self):
        """Set default values for the constrainedAstrometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "simpleFlux"
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        # See Readme for an explanation of these empirical values.
        pa1 = 0.024
        metrics = {'collected_photometry_refStars': 2948,
                   'selected_photometry_refStars': 657,
                   'associated_photometry_fittedStars': 1151,
                   'selected_photometry_fittedStars': 851,
                   'selected_photometry_ccdImages': 6,
                   'photometry_final_chi2': 1811.45,
                   'photometry_final_ndof': 1072
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_simpleFlux(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self.config.photometryModel = "simpleMagnitude"
        metrics['photometry_final_chi2'] = 1845.67
        metrics['photometry_final_ndof'] = 1074

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms(self):
        """Test that colorterms are applied and change the fit."""
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self.config.photometryModel = "simpleMagnitude"
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/hsc-colorterms-config.py')
        self.configfiles.append(test_config)

        # Final chi2 should be different, but I don't have an a-priori reason
        # to expect it to be larger or smaller.
        metrics['photometry_final_chi2'] = 1966.67

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms_no_library(self):
        """Fail Config validation if the color term library isn't defined."""
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/hsc-colorterms_no_library-config.py')
        self.configfiles.append(test_config)

        with self.assertRaises(lsst.pex.config.FieldValidationError):
            self._testJointcalTask(2, None, None, pa1)

    def testJointcalTask_2_visits_simple_astrometry_no_photometry(self):
        """Test turning off fitting photometry."""
        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2948,
                   'selected_astrometry_refStars': 657,
                   'associated_astrometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 851,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 819.07,
                   'astrometry_final_ndof': 2134,
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        data_refs = self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute,
                                           None, metrics=metrics)

        for data_ref in data_refs:
            with self.assertRaises(lsst.daf.persistence.butlerExceptions.NoResults):
                data_ref.get('jointcal_photoCalib')

    def test_jointcalTask_2_visits_simple_astrometry_gaia_refcat(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        # use the a.net refcat for photometry, gaia for astrometry
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-gaia-config.py')
        self.configfiles.append(test_config)
        dist_rms_relative = 17e-3*u.arcsecond

        # See Readme for an explanation of these empirical values.
        # NOTE: PA1 is slightly different here, because the number of SDSS
        # cross-matches within 0.1" goes down after we apply the GAIA-fit WCS.
        pa1 = 0.02405
        metrics = {'collected_astrometry_refStars': 1425,
                   'collected_photometry_refStars': 2948,
                   'selected_astrometry_refStars': 271,
                   'selected_photometry_refStars': 657,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 645,
                   'selected_photometry_fittedStars': 851,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_final_chi2': 435.01995,
                   'astrometry_final_ndof': 1412,
                   'photometry_final_chi2': 1811.45,
                   'photometry_final_ndof': 1072
                   }
        # NOTE: The astrometry/photometry tests are computed using the a.net SDSS refcat,
        # so the absolute astrometry RMS will be larger (because GAIA is better, so
        # comparing against SDSS will look "worse").
        dist_rms_absolute = 56e-3*u.arcsecond
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simple_astrometry_no_photometry_match_cut_10(self):
        """A larger matching radius will result in more associated fittedStars.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.matchCut = 10.0  # TODO: once DM-6885 is fixed, we need to put `*lsst.afw.geom.arcseconds`
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2948,
                   'selected_astrometry_refStars': 693,
                   'associated_astrometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 876,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 818.46,
                   'astrometry_final_ndof': 2132,
                   }
        pa1 = None
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry(self):
        """3 visit, default config to compare with min_measurements_3 test."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.minMeasurements = 2
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2948,
                   'selected_astrometry_refStars': 695,
                   'associated_astrometry_fittedStars': 1270,
                   'selected_astrometry_fittedStars': 1011,
                   'selected_astrometry_ccdImages': 8,
                   'astrometry_final_chi2': 1357.59,
                   'astrometry_final_ndof': 3302,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry_min_measurements_3(self):
        """Raising min_measurements to 3 will reduce the number of selected
        fitted stars (and thus the chisq and Ndof), but should not change the
        other values."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.minMeasurements = 3
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2948,
                   'selected_astrometry_refStars': 695,
                   'associated_astrometry_fittedStars': 1270,
                   'selected_astrometry_fittedStars': 808,
                   'selected_astrometry_ccdImages': 8,
                   'astrometry_final_chi2': 1210.05,
                   'astrometry_final_ndof': 2906,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
