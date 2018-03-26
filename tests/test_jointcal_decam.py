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


class JointcalTestDECAM(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'decam_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # See Readme for an explanation of this empirical value.
        self.dist_rms_absolute = 62.5e-3*u.arcsecond

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

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")

        # See Readme for an explanation of these empirical values.
        relative_error = 19e-3*u.arcsecond
        pa1 = 0.14
        metrics = {'collected_astrometry_refStars': 4866,
                   'collected_photometry_refStars': 4865,
                   'selected_astrometry_refStars': 661,
                   'selected_photometry_refStars': 661,
                   'associated_astrometry_fittedStars': 6749,
                   'associated_photometry_fittedStars': 6749,
                   'selected_astrometry_fittedStars': 2044,
                   'selected_photometry_fittedStars': 2044,
                   'selected_astrometry_ccdImages': 14,
                   'selected_photometry_ccdImages': 14,
                   'astrometry_final_chi2': 1974.13,
                   'astrometry_final_ndof': 3822,
                   'photometry_final_chi2': 3453.16,
                   'photometry_final_ndof': 2079,
                   }

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        """Help keep the two constrainedAstrometry tests consistent and make
        the difference between them more obvious."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryModel = "constrained"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # See Readme for an explanation of these empirical values.
        relative_error = 17e-3*u.arcsecond
        pa1 = None
        metrics = {'collected_astrometry_refStars': 4866,
                   'selected_astrometry_refStars': 661,
                   'associated_astrometry_fittedStars': 6749,
                   'selected_astrometry_fittedStars': 2044,
                   'selected_astrometry_ccdImages': 14,
                   'astrometry_final_chi2': 2072.89,
                   'astrometry_final_ndof': 3970,
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

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def setup_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        """Help keep the constrainedPhotometry tests consistent and make
        the differences between them more obvious.
        """
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.photometryModel = "constrained"
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

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
    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
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


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
