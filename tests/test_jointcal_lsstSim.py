# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from builtins import range

import inspect
import unittest
import os

from astropy import units as u

import lsst.afw.geom
import lsst.afw.image
import lsst.afw.coord
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase
import lsst.jointcal.jointcal

try:
    data_dir = lsst.utils.getPackageDir('testdata_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'twinkles1_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
dist_rms_absolute = 42e-3*u.arcsecond


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestLSSTSim(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    def setUp(self):
        do_plot = False

        # position of the Twinkles run 1 catalog
        center = lsst.afw.coord.IcrsCoord(53.00914*lsst.afw.geom.degrees, -27.43895*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(data_dir, 'twinkles1')
        all_visits = list(range(840, 850))
        other_args = ['raft=2,2', 'sensor=1,1', 'filter=r']

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    @unittest.skip('jointcal currently fails (may segfault) if only given one catalog!')
    def testJointcalTask_1_visits(self):
        dist_rms_relative = 0*u.arcsecond  # there is no such thing as a "relative" test for 1 catalog.
        pa1 = 2.64e-3
        self._testJointcalTask(1, dist_rms_relative, dist_rms_absolute, pa1)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_2_visits(self):
        # NOTE: The relative RMS limits were empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 9.7e-3*u.arcsecond
        pa1 = 2.64e-3
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_10_visits(self):
        dist_rms_relative = 7.9e-3*u.arcsecond
        pa1 = 2.64e-3
        self._testJointcalTask(10, dist_rms_relative, dist_rms_absolute, pa1)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_2_visits_no_astrometry(self):
        """Test turning off fitting astrometry."""
        pa1 = 2.64e-3
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller)
        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)

        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        self.assertIsNone(rms_result.dist_relative)
        self.assertIsNone(rms_result.dist_absolute)
        self.assertLess(rms_result.pa1, pa1)

        for data_ref in data_refs:
            wcs = data_ref.get('wcs').getWcs()
            self.assertIsNone(wcs)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_2_visits_no_photometry(self):
        """Test turning off fitting photometry."""
        dist_rms_relative = 9.7e-3*u.arcsecond
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller)
        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)

        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        self.assertLess(rms_result.dist_relative, dist_rms_relative)
        self.assertLess(rms_result.dist_absolute, dist_rms_absolute)
        self.assertIsNone(rms_result.pa1)

        for data_ref in data_refs:
            calib = data_ref.get('wcs').getCalib()
            blank_calib = lsst.afw.image.Calib()
            self.assertEqual(calib, blank_calib)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
