# See COPYRIGHT file at the top of the source tree.
import unittest
import os

from astropy import units as u

import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHT(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'cfht_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # We don't want the absolute astrometry to become significantly worse
        # than the single-epoch astrometry (about 0.040").
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 48.6e-3*u.arcsecond

        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.geom.SpherePoint(214.884832, 52.6622199, lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'cfht')
        all_visits = [849375, 850587]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot,
                        log_level="DEBUG")

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")

        # to test whether we got the expected chi2 contribution files.
        self.other_args.extend(['--config', 'writeChi2ContributionFiles=True'])

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 11e-3*u.arcsecond
        pa1 = 0.014
        metrics = {'collected_astrometry_refStars': 825,
                   'collected_photometry_refStars': 825,
                   'selected_astrometry_refStars': 350,
                   'selected_photometry_refStars': 350,
                   'associated_astrometry_fittedStars': 2269,
                   'associated_photometry_fittedStars': 2269,
                   'selected_astrometry_fittedStars': 1239,
                   'selected_photometry_fittedStars': 1239,
                   'selected_astrometry_ccdImages': 12,
                   'selected_photometry_ccdImages': 12,
                   'astrometry_final_chi2': 1150.62,
                   'astrometry_final_ndof': 2550,
                   'photometry_final_chi2': 2824.86,
                   'photometry_final_ndof': 1388
                   }

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

        # Check for the existence of the chi2 contribution files.
        expected = ['photometry_initial_chi2-0_r', 'astrometry_initial_chi2-0_r',
                    'photometry_final_chi2-0_r', 'astrometry_final_chi2-0_r']
        for partial in expected:
            name = partial+'-ref.csv'
            self.assertTrue(os.path.exists(name), msg="Did not find file %s"%name)
            os.remove(name)
            name = partial+'-meas.csv'
            self.assertTrue(os.path.exists(name), msg='Did not find file %s'%name)
            os.remove(name)

    def setup_jointcalTask_2_visits_constrainedAstrometry(self):
        """Help keep the constrainedAstrometry tests consistent and make
        the differences between them more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryModel = "constrained"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        dist_rms_relative = 12e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 825,
                   'selected_astrometry_refStars': 350,
                   'associated_astrometry_fittedStars': 2269,
                   'selected_astrometry_fittedStars': 1239,
                   'selected_astrometry_ccdImages': 12,
                   'astrometry_final_chi2': 1253.80,
                   'astrometry_final_ndof': 2660,
                   }

        return dist_rms_relative, metrics

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        dist_rms_relative, metrics = self.setup_jointcalTask_2_visits_constrainedAstrometry()
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedAstrometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't affect astrometry.
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
        metrics['astrometry_final_chi2'] = 922.76
        metrics['astrometry_final_ndof'] = 2486

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, None, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedPhotometry(self):
        """Help keep the constrainedPhotometry tests consistent and make
        the differences between them more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.photometryModel = "constrained"
        self.config.doAstrometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_astrometry = False

        # See Readme for an explanation of these empirical values.
        pa1 = 0.017
        metrics = {'collected_photometry_refStars': 825,
                   'selected_photometry_refStars': 350,
                   'associated_photometry_fittedStars': 2269,
                   'selected_photometry_fittedStars': 1239,
                   'selected_photometry_ccdImages': 12,
                   'photometry_final_chi2': 2655.86,
                   'photometry_final_ndof': 1328
                   }
        return pa1, metrics

    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_no_rank_update(self):
        """Demonstrate that skipping the rank update doesn't affect photometry.
        """
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.photometryDoRankUpdate = False

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_lineSearch(self):
        """Activating the line search should only slightly change the chi2 in this case."""
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.allowLineSearch = True

        # Only this value should differ from the metrics defined in setup above.
        metrics['photometry_final_chi2'] = 2642.47

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_flagged(self):
        """Test the use of the FlaggedSourceSelector."""
        pa1, metrics = self.setup_jointcalTask_2_visits_constrainedPhotometry()
        self.config.sourceSelector.name = "flagged"
        # Reduce warnings due to flaggedSourceSelector having fewer sources than astrometrySourceSelector.
        self.config.minMeasuredStarsPerCcd = 30
        self.config.minRefStarsPerCcd = 20

        pa1 = 0.026
        metrics['selected_photometry_refStars'] = 212
        metrics['associated_photometry_fittedStars'] = 270
        metrics['selected_photometry_fittedStars'] = 244
        metrics['photometry_final_chi2'] = 369.96
        metrics['photometry_final_ndof'] = 252

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
