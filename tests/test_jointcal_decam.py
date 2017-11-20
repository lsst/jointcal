# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from builtins import str
from builtins import range

import unittest
import os

from astropy import units as u

import lsst.afw.coord
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
        center = lsst.afw.coord.IcrsCoord(150.1191666*lsst.afw.geom.degrees, 2.20583333*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'decam')
        all_visits = [176837, 176846]
        ccdnums = '^'.join(str(x) for x in range(10, 19))
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

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 32e-3*u.arcsecond
        pa1 = 0.14
        # NOTE: decam fits are currently not converging; the chi2 jumps around, so skip that Metric.
        metrics = {'collected_astrometry_refStars': 8194,
                   'collected_photometry_refStars': 8190,
                   'selected_astrometry_refStars': 8194,
                   'selected_photometry_refStars': 8190,
                   'associated_astrometry_fittedStars': 8241,
                   'associated_photometry_fittedStars': 8241,
                   'selected_astrometry_fittedStars': 2261,
                   'selected_photometry_fittedStars': 2261,
                   'selected_astrometry_ccdImages': 17,
                   'selected_photometry_ccdImages': 17,
                   'astrometry_final_chi2': None,
                   'astrometry_final_ndof': 4306,
                   'photometry_final_chi2': None,
                   'photometry_final_ndof': 2333,
                   }

        # TODO DM-12653: decam fails due to factorization problems.
        with self.assertRaises(RuntimeError):
            self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPoly(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryModel = "constrainedPoly"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 20e-3*u.arcsecond
        pa1 = None
        metrics = {'collected_astrometry_refStars': 8194,
                   'selected_astrometry_refStars': 8194,
                   'associated_astrometry_fittedStars': 8241,
                   'selected_astrometry_fittedStars': 2261,
                   'selected_astrometry_ccdImages': 17,
                   'astrometry_final_chi2': 5106.2,
                   'astrometry_final_ndof': 4530,
                   }

        # TODO DM-12653: decam fails due to factorization problems.
        with self.assertRaises(RuntimeError):
            self._testJointcalTask(2, relative_error, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
