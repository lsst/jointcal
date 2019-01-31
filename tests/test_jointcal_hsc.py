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


class JointcalTestHSC(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'hsc_and_index')
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

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.024
        metrics = {'collected_astrometry_refStars': 2187,
                   'collected_photometry_refStars': 2187,
                   'selected_astrometry_refStars': 515,
                   'selected_photometry_refStars': 515,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 770,
                   'selected_photometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_final_chi2': 691.12,
                   'astrometry_final_ndof': 1858,
                   'photometry_final_chi2': 1557.27,
                   'photometry_final_ndof': 968
                   }
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_11_visits_no_photometry(self):
        """Test 11 visits with different filters.
        Testing photometry doesn't make sense for this currently.
        """

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = None  # pa1 = 0.134
        metrics = {'collected_astrometry_refStars': 3649,
                   'selected_astrometry_refStars': 1038,
                   'associated_astrometry_fittedStars': 2908,
                   'selected_astrometry_fittedStars': 2203,
                   'selected_astrometry_ccdImages': 33,
                   'astrometry_final_chi2': 7929.656,
                   'astrometry_final_ndof': 14262,
                   }
        self._testJointcalTask(11, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_simplePhotometry(self):
        """Help keep the simplePhotometry tests consistent and make
        the differences between them more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "simpleFlux"
        self.config.doAstrometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_astrometry = False

        # See Readme for an explanation of these empirical values.
        pa1 = 0.024
        metrics = {'collected_photometry_refStars': 2187,
                   'selected_photometry_refStars': 515,
                   'associated_photometry_fittedStars': 1151,
                   'selected_photometry_fittedStars': 770,
                   'selected_photometry_ccdImages': 6,
                   'photometry_final_chi2': 1557.27,
                   'photometry_final_ndof': 968
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_simpleFlux(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simpleMagnitude(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_simplePhotometry()
        self.config.photometryModel = "simpleMagnitude"
        metrics['photometry_final_chi2'] = 1550.47
        metrics['photometry_final_ndof'] = 967

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def testJointcalTask_2_visits_no_photometry(self):
        """Test turning off fitting photometry."""
        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2187,
                   'selected_astrometry_refStars': 515,
                   'associated_astrometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 691.1210,
                   'astrometry_final_ndof': 1858,
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        data_refs = self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute,
                                           None, metrics=metrics)

        for data_ref in data_refs:
            with self.assertRaises(lsst.daf.persistence.butlerExceptions.NoResults):
                data_ref.get('jointcal_photoCalib')

    def test_jointcalTask_2_visits_gaia_refcat(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        # use the a.net refcat for photometry, gaia for astrometry
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-config.py')
        self.configfiles.extend([test_config])
        dist_rms_relative = 17e-3*u.arcsecond

        # See Readme for an explanation of these empirical values.
        # NOTE: PA1 is slightly different here, because the number of SDSS
        # cross-matches within 0.1" goes down after we apply the GAIA-fit WCS.
        pa1 = 0.02405
        metrics = {'collected_astrometry_refStars': 1425,
                   'collected_photometry_refStars': 2187,
                   'selected_astrometry_refStars': 271,
                   'selected_photometry_refStars': 515,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 645,
                   'selected_photometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_final_chi2': 435.01995,
                   'astrometry_final_ndof': 1412,
                   'photometry_final_chi2': 1557.27,
                   'photometry_final_ndof': 968
                   }
        # NOTE: The astrometry/photometry tests are computed using the a.net SDSS refcat,
        # so the absolute astrometry RMS will be larger (because GAIA is better, so
        # comparing against SDSS will look "worse").
        dist_rms_absolute = 56e-3*u.arcsecond
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_no_photometry_match_cut_10(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.matchCut = 10.0  # TODO: once DM-6885 is fixed, we need to put `*lsst.afw.geom.arcseconds`
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2187,
                   'selected_astrometry_refStars': 546,
                   'associated_astrometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 790,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_final_chi2': 690.509,
                   'astrometry_final_ndof': 1856,
                   }
        pa1 = None
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_no_photometry(self):
        """3 visit, default config to compare with min_measurements_3 test."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.minMeasurements = 2
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2187,
                   'selected_astrometry_refStars': 541,
                   'associated_astrometry_fittedStars': 1270,
                   'selected_astrometry_fittedStars': 946,
                   'selected_astrometry_ccdImages': 8,
                   'astrometry_final_chi2': 1229.212,
                   'astrometry_final_ndof': 3008,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_no_photometry_min_measurements_3(self):
        """Raising min_measurements to 3 will reduce the number of selected
        fitted stars (and thus the chisq and Ndof), but should not change the
        other values."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.minMeasurements = 3
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2187,
                   'selected_astrometry_refStars': 541,
                   'associated_astrometry_fittedStars': 1270,
                   'selected_astrometry_fittedStars': 696,
                   'selected_astrometry_ccdImages': 8,
                   'astrometry_final_chi2': 1047.57,
                   'astrometry_final_ndof': 2526,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
