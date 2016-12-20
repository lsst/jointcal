# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os

from astropy import units as u

import lsst.afw.coord
import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('validation_data_hsc')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'sdss-dr9-fink-v5b')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
absolute_error = 52e-3*u.arcsecond


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestHSC(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    def setUp(self):
        do_plot = False

        # center of the hsc validation_data catalog
        center = lsst.afw.coord.IcrsCoord(320.367492*lsst.afw.geom.degrees, 0.3131554*lsst.afw.geom.degrees)
        radius = 5*lsst.afw.geom.degrees

        input_dir = os.path.join(data_dir, 'DATA', 'rerun', 'validate_drp')
        # self.all_visits = [903334, 903336, 903338, 903342, 903344, 903346]
        all_visits = [903982, 904006, 904828, 904846]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

    @unittest.skipIf(data_dir is None, "validation_data_hsc not setup")
    def test_jointcalTask_2_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 17e-3*u.arcsecond
        self._testJointCalTask(2, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_hsc not setup")
    def test_jointcalTask_4_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 17e-3*u.arcsecond
        self._testJointCalTask(4, relative_error, absolute_error)

    # @unittest.skipIf(data_dir is None, "validation_data_hsc not setup")
    # def test_jointcalTask_6_visits(self):
    #     # NOTE: The relative RMS limit was empirically determined from the
    #     # first run of jointcal on this data. We should always do better than
    #     # this in the future!
    #     relative_error = 10e-3*u.arcsecond
    #     self._testJointCalTask(6, relative_error, absolute_error)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
