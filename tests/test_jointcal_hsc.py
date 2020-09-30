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
        except LookupError:
            raise unittest.SkipTest("testdata_jointcal not setup")
        try:
            lsst.utils.getPackageDir('obs_subaru')
        except LookupError:
            raise unittest.SkipTest("obs_subaru not setup")

    def setUp(self):
        # See Readme for an explanation of these empirical values.
        self.dist_rms_absolute = 21e-3*u.arcsecond
        self.dist_rms_relative = 7e-3*u.arcsecond

        do_plot = False

        # center of the hsc validation_data catalog
        center = lsst.geom.SpherePoint(337.710899, +0.807006, lsst.geom.degrees)
        radius = 0.5*lsst.geom.degrees

        input_dir = os.path.join(self.data_dir, 'hsc')
        all_visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

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
        pa1 = 0.016
        metrics = {'collected_astrometry_refStars': 568,
                   'collected_photometry_refStars': 6485,
                   'selected_astrometry_refStars': 137,
                   'selected_photometry_refStars': 1609,
                   'associated_astrometry_fittedStars': 2070,
                   'associated_photometry_fittedStars': 2070,
                   'selected_astrometry_fittedStars': 989,
                   'selected_photometry_fittedStars': 1731,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_final_chi2': 683.783,
                   'astrometry_final_ndof': 1916,
                   'photometry_final_chi2': 4997.62,
                   'photometry_final_ndof': 2188
                   }
        self._testJointcalTask(2, self.dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_10_visits_simple_astrometry_no_photometry(self):
        """Test all 10 visits with different filters.
        Testing photometry doesn't make sense for this currently.
        """

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_absolute = 23e-3*u.arcsecond
        dist_rms_relative = 13e-3*u.arcsecond
        pa1 = None
        metrics = {'collected_astrometry_refStars': 1316,
                   'selected_astrometry_refStars': 318,
                   'associated_astrometry_fittedStars': 5860,
                   'selected_astrometry_fittedStars': 3568,
                   'selected_astrometry_ccdImages': 30,
                   'astrometry_final_chi2': 9621.29,
                   'astrometry_final_ndof': 18562,
                   }
        self._testJointcalTask(10, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_simplePhotometry(self):
        """Set default values for the simplePhotometry tests, and make it so
        the differences between each test and the defaults are more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "simpleFlux"
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        # See Readme for an explanation of these empirical values.
        pa1 = 0.016
        metrics = {'collected_photometry_refStars': 6485,
                   'selected_photometry_refStars': 1609,
                   'associated_photometry_fittedStars': 2070,
                   'selected_photometry_fittedStars': 1731,
                   'selected_photometry_ccdImages': 6,
                   'photometry_final_chi2': 4997.62,
                   'photometry_final_ndof': 2188
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_simpleFlux(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self.config.photometryModel = "simpleMagnitude"
        metrics['photometry_final_chi2'] = 5236.91
        metrics['photometry_final_ndof'] = 2200

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms(self):
        """Test that colorterms are applied and change the fit."""
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self.config.photometryModel = "simpleMagnitude"
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/hsc-colorterms-config.py')
        self.configfiles.append(test_config)

        # slightly fewer refstars because they each need the specific filters required by the colorterms
        metrics['collected_photometry_refStars'] = 6478
        # Final chi2 should be different, but I don't have an a-priori reason
        # to expect it to be larger or smaller.
        metrics['photometry_final_chi2'] = 5181.25
        metrics['photometry_final_ndof'] = 2197

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms_no_library(self):
        """Fail Config validation if the color term library isn't defined."""
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                   'tests/config/hsc-colorterms_no_library-config.py')
        self.configfiles.append(test_config)

        with self.assertRaises(lsst.pex.config.FieldValidationError):
            self._testJointcalTask(2, None, None, pa1)

    def test_JointcalTask_2_visits_simple_astrometry_no_photometry(self):
        """Test turning off fitting photometry."""
        metrics = {'collected_astrometry_refStars': 568,
                   'selected_astrometry_refStars': 137,
                   'associated_astrometry_fittedStars': 2070,
                   'selected_astrometry_fittedStars': 989,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 683.783,
                   'astrometry_final_ndof': 1916,
                   }
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        data_refs = self._testJointcalTask(2, self.dist_rms_relative, self.dist_rms_absolute,
                                           None, metrics=metrics)

        for data_ref in data_refs:
            with self.assertRaises(lsst.daf.persistence.butlerExceptions.NoResults):
                data_ref.get('jointcal_photoCalib')

    def test_jointcalTask_2_visits_simple_astrometry_no_photometry_match_cut_10(self):
        """A larger matching radius will result in more associated fittedStars,
        and a somewhat worse final fit because stars that should not have been
        associated, were.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.matchCut = 10.0  # TODO: once DM-6885 is fixed, we need to put `*lsst.geom.arcseconds`
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # Slightly larger absolute astrometry RMS because of the larger matching radius
        dist_rms_absolute = 23e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 568,
                   'selected_astrometry_refStars': 211,
                   'associated_astrometry_fittedStars': 2070,
                   'selected_astrometry_fittedStars': 1042,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 683.369,
                   'astrometry_final_ndof': 1910,
                   }
        pa1 = None
        self._testJointcalTask(2, self.dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry(self):
        """3 visit, default config to compare with min_measurements_3 test."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.minMeasurements = 2
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # More visits and slightly worse relative and absolute rms.
        dist_rms_relative = 8.2e-3*u.arcsecond
        dist_rms_absolute = 24e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 1038,
                   'selected_astrometry_refStars': 209,
                   'associated_astrometry_fittedStars': 3199,
                   'selected_astrometry_fittedStars': 1282,
                   'selected_astrometry_ccdImages': 9,
                   'astrometry_final_chi2': 1013.20,
                   'astrometry_final_ndof': 2900,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry_min_measurements_3(self):
        """Raising min_measurements to 3 will reduce the number of selected
        fitted stars (and thus the chisq and Ndof), but should not change the
        other values."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.minMeasurements = 3
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # More visits and slightly worse relative and absolute rms.
        dist_rms_relative = 11e-3*u.arcsecond
        dist_rms_absolute = 24e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 1038,
                   'selected_astrometry_refStars': 209,
                   'associated_astrometry_fittedStars': 3199,
                   'selected_astrometry_fittedStars': 432,
                   'selected_astrometry_ccdImages': 9,
                   'astrometry_final_chi2': 360.884,
                   'astrometry_final_ndof': 1294,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
