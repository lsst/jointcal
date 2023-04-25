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
        input_dir = os.path.join(self.data_dir, 'cfht')
        all_visits = [849375, 850587]

        where = "instrument='MegaPrime' and tract=0 and skymap='discrete'"
        inputCollections = ["singleFrame", "skymaps"]
        refcats = {"gaia_dr2_20200414": os.path.join(input_dir, "gaia_dr2_20200414.ecsv"),
                   "ps1_pv3_3pi_20170110": os.path.join(input_dir, "ps1_pv3_3pi_20170110.ecsv"),
                   "sdss_dr9_fink_v5b": os.path.join(input_dir, "sdss-dr9-fink-v5b.ecsv")}

        outputDataId = {'instrument': 'MegaPrime', 'tract': 0, 'skymap': 'discrete'}
        self.setUp_base("lsst.obs.cfht.MegaPrime", "MegaPrime",
                        input_dir=input_dir,
                        all_visits=all_visits,
                        where=where,
                        inputCollections=inputCollections,
                        refcats=refcats,
                        refcatPath=input_dir,
                        outputDataId=outputDataId,
                        log_level="DEBUG")

        # The CFHT tests all produce the same set of output visits+detectors,
        # whether astrometry or photometry.
        self.outputVisits = {849375: (12, 13, 14, 21, 22, 23),
                             850587: (12, 13, 14, 21, 22, 23)}

    def test_jointcalTask_2_visits_simple(self):
        """Test the simple models with two visits and check that some debug
        output files also get created.
        """
        configOptions = {"astrometryModel": "simple", "photometryModel": "simpleFlux",
                         "writeInitialModel": True, "writeChi2FilesInitialFinal": True}

        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            configOptions["debugOutputPath"] = tempdir
            metrics = {'astrometry_collected_refStars': 867,
                       'photometry_collected_refStars': 11569,
                       'astrometry_prepared_refStars': 323,
                       'photometry_prepared_refStars': 2302,
                       'astrometry_matched_fittedStars': 2399,
                       'photometry_matched_fittedStars': 2399,
                       'astrometry_prepared_fittedStars': 1255,
                       'photometry_prepared_fittedStars': 2317,
                       'astrometry_prepared_ccdImages': 12,
                       'photometry_prepared_ccdImages': 12,
                       'astrometry_final_chi2': 1568.76,
                       'astrometry_final_ndof': 2356,
                       'photometry_final_chi2': 11561.1,
                       'photometry_final_ndof': 2849
                       }
            self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                                  astrometryOutputs=self.outputVisits, photometryOutputs=self.outputVisits)

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
        configOptions = {"astrometryModel": "constrained", "doPhotometry": False}
        metrics = {'astrometry_collected_refStars': 867,
                   'astrometry_prepared_refStars': 323,
                   'astrometry_matched_fittedStars': 2399,
                   'astrometry_prepared_fittedStars': 1255,
                   'astrometry_prepared_ccdImages': 12,
                   'astrometry_final_chi2': 1611.57,
                   'astrometry_final_ndof': 2428,
                   }
        return configOptions, metrics

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        configOptions['writeInitialModel'] = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            configOptions['debugOutputPath'] = tempdir

            self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                                  astrometryOutputs=self.outputVisits)

            filename = os.path.join(tempdir, "initial_astrometry_model-0_r.MP9601.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedAstrometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect astrometry.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        metrics['astrometry_final_chi2'] = 1611.57
        metrics['astrometry_final_ndof'] = 2318

        configOptions['astrometryDoRankUpdate'] = False

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedAstrometry_4sigma_outliers(self):
        """4 sigma outlier rejection means fewer available sources after the
        fitter converges, resulting in a smaller ndof and chi2.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        configOptions['outlierRejectSigma'] = 4
        metrics['astrometry_final_chi2'] = 1173.75
        metrics['astrometry_final_ndof'] = 2254

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryOutlierRelativeTolerance(self):
        """Test that astrometryOutlierRelativeTolerance changes the fit. Setting
        1% for the astrometryOutlierRelativeTolerance will result in higher chi2
        and ndof.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        configOptions['astrometryOutlierRelativeTolerance'] = 0.01
        metrics['astrometry_final_chi2'] = 2229.21
        metrics['astrometry_final_ndof'] = 2552

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryReferenceUncertainty_smaller(self):
        """Test with a smaller fake reference uncertainty: chi2 will be higher."""
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        astrometryRefErrConfig = os.path.join(self.path, 'config/astrometryReferenceErr-config.py')
        metrics['astrometry_final_chi2'] = 1479.02
        metrics['astrometry_final_ndof'] = 2524

        self._runJointcalTest(configFiles=[astrometryRefErrConfig],
                              configOptions=configOptions, metrics=metrics,
                              astrometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedAstrometry_astrometryReferenceUncertainty_None_fails(self):
        """Setting astrometryReferenceUncertainty=None should fail for refcats
        that have no position errors.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        badRefErrConfig = os.path.join(self.path, 'config/astrometryReferenceErr-None-config.py')
        with self.assertRaisesRegex(lsst.pex.config.FieldValidationError,
                                    "Reference catalog does not contain coordinate errors"):
            self._runJointcalTest(configFiles=[badRefErrConfig], configOptions=configOptions, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedPhotometry(self):
        """Set default values for the constrainedPhotometry tests, and make
        the differences between each test and the defaults more obvious.
        """
        configOptions = {"photometryModel": "constrainedFlux", "doAstrometry": False}

        metrics = {'photometry_collected_refStars': 11569,
                   'photometry_prepared_refStars': 2302,
                   'photometry_matched_fittedStars': 2399,
                   'photometry_prepared_fittedStars': 2317,
                   'photometry_prepared_ccdImages': 12,
                   'photometry_final_chi2': 11264.28,
                   'photometry_final_ndof': 2821
                   }
        return configOptions, metrics

    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        configOptions['writeInitialModel'] = True  # write the initial models
        # use a temporary directory for debug output, to prevent test collisions
        with tempfile.TemporaryDirectory() as tempdir:
            configOptions['debugOutputPath'] = tempdir

            self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                                  photometryOutputs=self.outputVisits)
            filename = os.path.join(tempdir, "initial_photometry_model-0_r.MP9601.txt")
            self.assertTrue(os.path.exists(filename), msg=f"Did not find file {filename}")

    def test_jointcalTask_2_visits_constrainedPhotometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't substantially affect photometry.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        configOptions['photometryDoRankUpdate'] = False

        # The constrainedPhotometry model is not purely linear, so a small
        # change in final chi2 is possible.
        metrics['photometry_final_chi2'] = 10896.76
        metrics['photometry_final_ndof'] = 2787

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              photometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedPhotometry_lineSearch(self):
        """Activating the line search should only slightly change the chi2.

        Activating line search for constrainedPhotometry should result in
        nearly the same final fit (the system is somewhat non-linear, so it
        may not be exactly the same: check the "Line search scale factor"
        lines in the DEBUG log for values that are not ~1 for proof).
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        configOptions['allowLineSearch'] = True

        metrics['photometry_final_chi2'] = 10000.87
        metrics['photometry_final_ndof'] = 2773

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              photometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedMagnitude_no_astrometry(self):
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        configOptions['photometryModel'] = "constrainedMagnitude"

        # The resulting fit should be close to the constrainedFlux model:
        # there are few CCDs and 2 visits, so there's not a lot of complexity
        # in this case to distinguish the flux vs. magnitude models.
        metrics['photometry_final_chi2'] = 10276.57
        metrics['photometry_final_ndof'] = 2823

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              photometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedFlux_pedestal(self):
        """Test that forcing a systematic flux error results in a lower chi2.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        # median fluxErr/flux in ps1 is 0.11, so we have to make this bigger
        # than that to actually allow more slop in the fit.
        configOptions['photometryErrorPedestal'] = 0.2

        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 3355.96
        # ndof may change; slightly different likelihood contours, and fewer
        # reference sources rejected.
        metrics['photometry_final_ndof'] = 3262

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              photometryOutputs=self.outputVisits)

    def test_jointcalTask_2_visits_constrainedMagnitude_pedestal(self):
        """Test that forcing a systematic flux error results in a lower chi2.
        """
        configOptions, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        configOptions['photometryModel'] = "constrainedMagnitude"
        # median fluxErr/flux in ps1 is 0.11, so we have to make this bigger
        # than that to actually allow more slop in the fit.
        configOptions['photometryErrorPedestal'] = 0.2

        # Final chi2 is much lower, because all sources contribute more error.
        metrics['photometry_final_chi2'] = 3165.87
        # ndof may change; slightly different likelihood contours, and fewer
        # reference sources rejected.
        metrics['photometry_final_ndof'] = 3224

        self._runJointcalTest(configOptions=configOptions, metrics=metrics,
                              photometryOutputs=self.outputVisits)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
