# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from builtins import range

import unittest
import os

from astropy import units as u

import lsst.afw.geom
import lsst.afw.coord
import lsst.utils
import lsst.pex.exceptions
import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('testdata_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'twinkles1_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
absolute_error = 42e-3*u.arcsecond


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
    @unittest.skip('jointcal currently fails if only given one catalog!')
    def testJointCalTask_1_visits(self):
        self._testJointCalTask(2, 0, absolute_error)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointCalTask_2_visits(self):
        # NOTE: The relative RMS limits were empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 9.7e-3*u.arcsecond
        self._testJointCalTask(2, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_4_visits(self):
        relative_error = 8.2e-3*u.arcsecond
        self._testJointCalTask(4, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_7_visits(self):
        relative_error = 8.1e-3*u.arcsecond
        self._testJointCalTask(7, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointCalTask_10_visits(self):
        relative_error = 7.9e-3*u.arcsecond
        self._testJointCalTask(10, relative_error, absolute_error)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
