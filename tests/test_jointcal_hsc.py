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
        input_dir = os.path.join(self.data_dir, 'hsc')
        all_visits = [34648, 34690, 34714, 34674, 34670, 36140, 35892, 36192, 36260, 36236]

        where = "instrument='HSC' and tract=9697 and skymap='hsc_rings_v1'"
        inputCollections = ["refcats/gen2",
                            "HSC/testdata",
                            "HSC/calib/unbounded"]
        outputDataId = {'instrument': 'HSC', 'tract': 9697, 'skymap': 'hsc_rings_v1'}

        self.setUp_base("lsst.obs.subaru.HyperSuprimeCam", "HSC",
                        input_dir=input_dir,
                        all_visits=all_visits,
                        where=where,
                        inputCollections=inputCollections,
                        outputDataId=outputDataId)

        self.configfiles.append(os.path.join(self.path, 'config/hsc-config.py'))

    def test_jointcalTask_2_visits_simple(self):
        configOptions = {"astrometryModel": "simple", "photometryModel": "simpleFlux"}
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        # test colorterm loading in gen3 (see DM-29884)
        configFiles = [os.path.join(self.path, 'config/hsc-colorterms-config.py')]

        # See Readme for an explanation of these empirical values.
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
        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64]}
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, configFiles=configFiles, metrics=metrics,
                              astrometryOutputs=outputVisits, photometryOutputs=outputVisits)

    def test_jointcalTask_ri_visits_2_bands_simple_dm32207(self):
        """Test gen3 butler jointcal putting 2 bands in same repo.
        Not testing metric values here; the goal of this test is to check that
        the metrics Connections are defined correctly so that the QauntumGraph
        generation doesn't fail.
        """
        configOptions = {"astrometryModel": "simple", "photometryModel": "simpleFlux"}
        where = (f" and visit in ({self.all_visits[0]},{self.all_visits[1]},"
                 f"{self.all_visits[5]},{self.all_visits[6]})")
        configFiles = []

        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64],
                        36140: [51, 59, 67], 35892: [23, 31, 39]}
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, configFiles=configFiles, metrics=None, nJobs=2,
                              astrometryOutputs=outputVisits, photometryOutputs=outputVisits)

    def test_jointcalTask_10_visits_simple_astrometry_no_photometry(self):
        """Test all 10 visits with different filters.
        Testing photometry doesn't make sense for this currently.
        """
        configOptions = {"astrometryModel": "simple", "doPhotometry": False}
        where = f" and visit in {tuple(self.all_visits)}"

        # See Readme for an explanation of these empirical values.
        # Gen3 jointcal splits the input by filter, so we get 15 ccdImages for
        # each of the r/i datasets (5 visits each). The metrics test only
        # checks the first Job (after sorting on the job keys), which is HSC-I.
        metrics = {'astrometry_collected_refStars': 897,
                   'astrometry_prepared_refStars': 192,
                   'astrometry_matched_fittedStars': 4238,
                   'astrometry_prepared_fittedStars': 2194,
                   'astrometry_prepared_ccdImages': 15,
                   'astrometry_final_chi2': 3557.92,
                   'astrometry_final_ndof': 7610,
                   }
        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64],
                        34714: [13, 18, 19], 34674: [15, 21, 28],
                        34670: [92, 93, 97], 36140: [51, 59, 67],
                        35892: [23, 31, 39], 36192: [8, 14, 20],
                        36260: [2, 7, 13], 36236: [87, 93, 98]}
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics, nJobs=2,
                              astrometryOutputs=outputVisits)

    def setup_jointcalTask_2_visits_simple_photometry(self):
        """Set default values for the simple_photometry tests, and make it so
        the differences between each test and the defaults are more obvious.
        """
        configOptions = {"photometryModel": "simpleFlux", "doAstrometry": False}
        # See Readme for an explanation of these empirical values.
        metrics = {'photometry_collected_refStars': 6485,
                   'photometry_prepared_refStars': 1609,
                   'photometry_matched_fittedStars': 2070,
                   'photometry_prepared_fittedStars': 1731,
                   'photometry_prepared_ccdImages': 6,
                   'photometry_final_chi2': 4997.62,
                   'photometry_final_ndof': 2188
                   }
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64]}
        return configOptions, metrics, where, outputVisits

    def test_jointcalTask_2_visits_simpleFlux(self):
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_photometry()
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              photometryOutputs=outputVisits)

    def test_jointcalTask_2_visits_simpleMagnitude(self):
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_photometry()
        configOptions['photometryModel'] = "simpleMagnitude"
        metrics['photometry_final_chi2'] = 5236.91
        metrics['photometry_final_ndof'] = 2200
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              photometryOutputs=outputVisits)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms(self):
        """Test that colorterms are applied and change the fit."""
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_photometry()
        configOptions['photometryModel'] = "simpleMagnitude"
        test_config = os.path.join(self.path, 'config/hsc-colorterms-config.py')

        # slightly fewer refstars because they each need the specific filters required by the colorterms
        metrics['photometry_collected_refStars'] = 6478
        # Final chi2 should be different, but I don't have an a-priori reason
        # to expect it to be larger or smaller.
        metrics['photometry_final_chi2'] = 5181.25
        metrics['photometry_final_ndof'] = 2197
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, configFiles=[test_config], metrics=metrics,
                              photometryOutputs=outputVisits)

    def test_jointcalTask_2_visits_simpleMagnitude_colorterms_no_library(self):
        """Fail Config validation if the color term library isn't defined."""
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_photometry()
        test_config = os.path.join(self.path, 'config/hsc-colorterms_no_library-config.py')
        self.configfiles.append(test_config)

        with self.assertRaises(lsst.pex.config.FieldValidationError):
            self._runJointcalTest(whereSuffix=where,
                                  configOptions=configOptions, configFiles=[test_config])

    def setup_jointcalTask_2_visits_simple_astrometry(self):
        """Set default values for the 2 visit simple_astrometry tests, so that
        the differences between each test and the defaults are more obvious.
        """
        configOptions = {"astrometryModel": "simple", "doPhotometry": False}
        # See Readme for an explanation of these empirical values.
        metrics = {'astrometry_collected_refStars': 568,
                   'astrometry_prepared_refStars': 137,
                   'astrometry_matched_fittedStars': 2070,
                   'astrometry_prepared_fittedStars': 989,
                   'astrometry_prepared_ccdImages': 6,
                   'astrometry_final_chi2': 835.473,
                   'astrometry_final_ndof': 1918,
                   }
        where = f" and visit in ({self.all_visits[0]},{self.all_visits[1]})"
        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64]}
        return configOptions, metrics, where, outputVisits

    def test_JointcalTask_2_visits_simple_astrometry_no_photometry(self):
        """Test turning off fitting photometry."""
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_astrometry()
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=outputVisits)

    def test_jointcalTask_2_visits_simple_astrometry_no_photometry_match_cut_10(self):
        """A larger matching radius will result in more associated fittedStars,
        and a somewhat worse final fit because stars that should not have been
        associated, were.
        """
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_2_visits_simple_astrometry()
        # TODO: once DM-6885 is fixed, we need to put `*lsst.geom.arcseconds`
        configOptions['matchCut'] = 10.0

        metrics['astrometry_prepared_refStars'] = 211
        metrics['astrometry_prepared_fittedStars'] = 1042
        metrics['astrometry_final_chi2'] = 819.608
        metrics['astrometry_final_ndof'] = 1872
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=outputVisits)

    def setup_jointcalTask_3_visits_simple_astrometry(self):
        """Set default values for the 3 visit simple_astrometry tests, so that
        the differences between each test and the defaults are more obvious.
        """
        configOptions = {"astrometryModel": "simple", "doPhotometry": False}
        # See Readme for an explanation of these empirical values.
        metrics = {'astrometry_collected_refStars': 1038,
                   'astrometry_prepared_refStars': 209,
                   'astrometry_matched_fittedStars': 3199,
                   'astrometry_prepared_fittedStars': 1282,
                   'astrometry_prepared_ccdImages': 9,
                   'astrometry_final_chi2': 1221.646,
                   'astrometry_final_ndof': 2892,
                   }
        where = f" and visit in {tuple(self.all_visits[:3])}"
        outputVisits = {34648: [51, 59, 67], 34690: [48, 56, 64], 34714: [13, 18, 19]}
        return configOptions, metrics, where, outputVisits

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry(self):
        """3 visit, default config to compare with min_measurements_3 test."""
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_3_visits_simple_astrometry()
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=outputVisits)

    def test_jointcalTask_3_visits_simple_astrometry_no_photometry_min_measurements_3(self):
        """Raising min_measurements to 3 will reduce the number of selected
        fitted stars (and thus the chisq and Ndof), but should not change the
        other values."""
        configOptions, metrics, where, outputVisits = self.setup_jointcalTask_3_visits_simple_astrometry()
        configOptions['minMeasurements'] = 3
        metrics['astrometry_prepared_fittedStars'] = 432
        metrics['astrometry_final_chi2'] = 567.451
        metrics['astrometry_final_ndof'] = 1286
        self._runJointcalTest(whereSuffix=where,
                              configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=outputVisits)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
