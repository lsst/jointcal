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
import numpy as np

from lsst.daf.butler import Butler
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

        where = "instrument='HSC' and tract=9697 and skymap='hsc_rings_v1'"
        inputCollections = ["refcats/gen2",
                            "HSC/testdata",
                            "HSC/calib/unbounded"]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot,
                        where=where,
                        inputCollections=inputCollections)

        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-config.py')
        self.configfiles.append(test_config)

    def test_jointcalTask_2_visits_simple(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.photometryModel = "simpleFlux"
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # See Readme for an explanation of these empirical values.
        pa1 = 0.016
        metrics = {'astrometry_collected_refStars': 568,
                   'photometry_collected_refStars': 6485,
                   'astrometry_prepared_refStars': 137,
                   'photometry_prepared_refStars': 1609,
                   'astrometry_matched_fittedStars': 2070,
                   'photometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 989,
                   'photometry_prepared_fittedStars': 1731,
                   'astrometry_prepared_ccdImages': 6,
                   'photometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 835.473,
                   'astrometry_final_ndof': 1918,
                   'photometry_final_chi2': 4997.62,
                   'photometry_final_ndof': 2188
                   }
        self._testJointcalTask(2, self.dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_simple_gen3(self):
        """Test gen3 butler jointcal."""
        configOptions = {"astrometryModel": "simple", "photometryModel": "simpleFlux"}
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        # test colorterm loading in gen3 (see DM-29884)
        colorterm_config = os.path.join(lsst.utils.getPackageDir('jointcal'),
                                        'tests/config/hsc-colorterms-config.py')
        configFiles = [os.path.join(lsst.utils.getPackageDir('jointcal'),
                                    'tests/config/config-gen3-hsc.py'),
                       colorterm_config]

        metrics = {'astrometry_collected_refStars': 568,
                   'photometry_collected_refStars': 6478,
                   'astrometry_prepared_refStars': 137,
                   'photometry_prepared_refStars': 1609,
                   'astrometry_matched_fittedStars': 2070,
                   'photometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 989,
                   'photometry_prepared_fittedStars': 1731,
                   'astrometry_prepared_ccdImages': 6,
                   'photometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 835.473,
                   'astrometry_final_ndof': 1918,
                   'photometry_final_chi2': 4977.2,
                   'photometry_final_ndof': 2188
                   }
        repo = self._runGen3Jointcal("lsst.obs.subaru.HyperSuprimeCam", "HSC", whereSuffix=where,
                                     configOptions=configOptions, configFiles=configFiles, metrics=metrics)
        # TODO DM-28863: this does not currently test anything other than the code
        # running without raising and that it writes non-empty output.
        butler = Butler(repo, collections=['HSC/tests/all'])

        def check_output(visit, detectors):
            """Check that there is something for each detector, and only
            entries for the correct detectors."""
            dataId = {'visit': visit, 'instrument': 'HSC', 'tract': 9697, 'skymap': 'hsc_rings_v1'}

            catalog = butler.get('jointcalPhotoCalibCatalog', dataId)
            for record in catalog:
                self.assertIsInstance(record.getPhotoCalib(), lsst.afw.image.PhotoCalib,
                                      msg=f"visit {visit}: {record}")
            np.testing.assert_array_equal(catalog['id'], detectors)

            catalog = butler.get('jointcalSkyWcsCatalog', dataId)
            for record in catalog:
                self.assertIsInstance(record.getWcs(), lsst.afw.geom.SkyWcs,
                                      msg=f"visit {visit}: {record}")
            np.testing.assert_array_equal(catalog['id'], detectors)

        check_output(34648, [51, 59, 67])
        check_output(34690, [48, 56, 64])

    def test_jointcalTask_ri_visits_2_bands_simple_gen3_dm32207(self):
        """Test gen3 butler jointcal putting 2 bands in same repo."""
        configOptions = {"astrometryModel": "simple", "photometryModel": "simpleFlux"}
        where = (f" and visit in ({self.all_visits[0]},{self.all_visits[1]},"
                 f"{self.all_visits[5]},{self.all_visits[6]})")
        configFiles = [os.path.join(lsst.utils.getPackageDir('jointcal'),
                                    'tests/config/config-gen3-hsc.py')]

        self._runGen3Jointcal("lsst.obs.subaru.HyperSuprimeCam", "HSC", whereSuffix=where,
                              configOptions=configOptions, configFiles=configFiles, metrics=None, nJobs=2)

    def test_jointcalTask_2_visits_simple_astrometry_no_photometry_gen3(self):
        """Test gen3 butler jointcal, no photometry."""
        configOptions = {"astrometryModel": "simple", "doPhotometry": False}
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        configFiles = [os.path.join(lsst.utils.getPackageDir('jointcal'),
                                    'tests/config/config-gen3-hsc.py')]

        metrics = {'astrometry_collected_refStars': 568,
                   'astrometry_prepared_refStars': 137,
                   'astrometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 989,
                   'astrometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 835.473,
                   'astrometry_final_ndof': 1918,
                   }
        repo = self._runGen3Jointcal("lsst.obs.subaru.HyperSuprimeCam", "HSC", whereSuffix=where,
                                     configOptions=configOptions, configFiles=configFiles, metrics=metrics)
        # TODO DM-28863: this does not currently test anything other than the code
        # running without raising and that it writes non-empty output.
        butler = Butler(repo, collections=['HSC/tests/all'])

        def check_output(visit, detectors):
            """Check that there is something for each detector, and only
            entries for the correct detectors."""
            dataId = {'visit': visit, 'instrument': 'HSC', 'tract': 9697, 'skymap': 'hsc_rings_v1'}

            catalog = butler.get('jointcalSkyWcsCatalog', dataId)
            for record in catalog:
                self.assertIsInstance(record.getWcs(), lsst.afw.geom.SkyWcs,
                                      msg=f"visit {visit}: {record}")
            np.testing.assert_array_equal(catalog['id'], detectors)

        check_output(34648, [51, 59, 67])
        check_output(34690, [48, 56, 64])

    def test_jointcalTask_2_visits_simple_photometry_no_astrometry_gen3(self):
        """Test gen3 butler jointcal, no astrometry."""
        configOptions = {"doAstrometry": False, "photometryModel": "simpleFlux"}
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        configFiles = [os.path.join(lsst.utils.getPackageDir('jointcal'),
                                    'tests/config/config-gen3-hsc.py')]

        metrics = {'photometry_collected_refStars': 6485,
                   'photometry_prepared_refStars': 1609,
                   'photometry_matched_fittedStars': 2070,
                   'photometry_prepared_fittedStars': 1731,
                   'photometry_prepared_ccdImages': 6,
                   'photometry_final_chi2': 4997.62,
                   'photometry_final_ndof': 2188
                   }
        repo = self._runGen3Jointcal("lsst.obs.subaru.HyperSuprimeCam", "HSC", whereSuffix=where,
                                     configOptions=configOptions, configFiles=configFiles, metrics=metrics)
        # TODO DM-28863: this does not currently test anything other than the code
        # running without raising and that it writes non-empty output.
        butler = Butler(repo, collections=['HSC/tests/all'])

        def check_output(visit, detectors):
            """Check that there is something for each detector, and only
            entries for the correct detectors."""
            dataId = {'visit': visit, 'instrument': 'HSC', 'tract': 9697, 'skymap': 'hsc_rings_v1'}

            catalog = butler.get('jointcalPhotoCalibCatalog', dataId)
            for record in catalog:
                self.assertIsInstance(record.getPhotoCalib(), lsst.afw.image.PhotoCalib,
                                      msg=f"visit {visit}: {record}")
            np.testing.assert_array_equal(catalog['id'], detectors)

        check_output(34648, [51, 59, 67])
        check_output(34690, [48, 56, 64])

    def test_jointcalTask_10_visits_simple_astrometry_no_photometry(self):
        """Test all 10 visits with different filters.
        Testing photometry doesn't make sense for this currently.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # See Readme for an explanation of these empirical values.
        dist_rms_absolute = 23e-3*u.arcsecond
        dist_rms_relative = 13e-3*u.arcsecond
        pa1 = None
        metrics = {'astrometry_collected_refStars': 1316,
                   'astrometry_prepared_refStars': 318,
                   'astrometry_matched_fittedStars': 5860,
                   'astrometry_prepared_fittedStars': 3568,
                   'astrometry_prepared_ccdImages': 30,
                   'astrometry_final_chi2': 10225.31,
                   'astrometry_final_ndof': 18576,
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
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # See Readme for an explanation of these empirical values.
        pa1 = 0.016
        metrics = {'photometry_collected_refStars': 6485,
                   'photometry_prepared_refStars': 1609,
                   'photometry_matched_fittedStars': 2070,
                   'photometry_prepared_fittedStars': 1731,
                   'photometry_prepared_ccdImages': 6,
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
        metrics['photometry_collected_refStars'] = 6478
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
        metrics = {'astrometry_collected_refStars': 568,
                   'astrometry_prepared_refStars': 137,
                   'astrometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 989,
                   'astrometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 835.473,
                   'astrometry_final_ndof': 1918,
                   }
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryModel = "simple"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

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
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # Slightly larger absolute astrometry RMS because of the larger matching radius
        dist_rms_absolute = 23e-3*u.arcsecond
        metrics = {'astrometry_collected_refStars': 568,
                   'astrometry_prepared_refStars': 211,
                   'astrometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 1042,
                   'astrometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 819.608,
                   'astrometry_final_ndof': 1872,
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
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # More visits and slightly worse relative and absolute rms.
        dist_rms_relative = 8.3e-3*u.arcsecond
        dist_rms_absolute = 24e-3*u.arcsecond
        metrics = {'astrometry_collected_refStars': 1038,
                   'astrometry_prepared_refStars': 209,
                   'astrometry_matched_fittedStars': 3199,
                   'astrometry_prepared_fittedStars': 1282,
                   'astrometry_prepared_ccdImages': 9,
                   'astrometry_final_chi2': 1221.63,
                   'astrometry_final_ndof': 2892,
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
        self.input_dir = os.path.join(self.data_dir, 'hsc/repo')

        # More visits and slightly worse relative and absolute rms.
        dist_rms_relative = 11e-3*u.arcsecond
        dist_rms_absolute = 24e-3*u.arcsecond
        metrics = {'astrometry_collected_refStars': 1038,
                   'astrometry_prepared_refStars': 209,
                   'astrometry_matched_fittedStars': 3199,
                   'astrometry_prepared_fittedStars': 432,
                   'astrometry_prepared_ccdImages': 9,
                   'astrometry_final_chi2': 567.433,
                   'astrometry_final_ndof': 1286,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
