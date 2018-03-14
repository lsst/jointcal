# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from builtins import str

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
        # This value was empirically determined from the first run of jointcal on
        # this data, and will likely vary from survey to survey.
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
                        do_plot=do_plot)

    def test_jointcalTask_2_visits(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 19e-3*u.arcsecond
        pa1 = 0.14
        # NOTE: decam fits are currently not converging; the chi2 jumps around, so skip that Metric.
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

    def test_jointcalTask_2_visits_constrainedAstrometry_no_photometry(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryModel = "constrained"
        self.config.doPhotometry = False
        self.config.sourceSelector['astrometry'].badFlags.append("base_PixelFlags_flag_interpolated")
        self.jointcalStatistics.do_photometry = False

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
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

        self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
