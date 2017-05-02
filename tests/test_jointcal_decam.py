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
        self.absolute_error = 62e-3*u.arcsecond

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
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 32e-3*u.arcsecond
        pa1 = 0.0383
        # NOTE: decam fits are currently not converging; the chi2 jumps around, so skip that test.
        metrics = {'collectedAstrometryRefStars': 8194,
                   'collectedPhotometryRefStars': 8194,
                   'selectedAstrometryRefStars': 8194,
                   'selectedPhotometryRefStars': 8194,
                   'associatedAstrometryFittedStars': 8241,
                   'associatedPhotometryFittedStars': 8241,
                   'selectedAstrometryFittedStars': 2261,
                   'selectedPhotometryFittedStars': 2261,
                   'selectedAstrometryCcdImageList': 17,
                   'selectedPhotometryCcdImageList': 17,
                   'astrometryFinalChi2': None,
                   'astrometryFinalNdof': 4306,
                   'photometryFinalChi2': None,
                   'photometryFinalNdof': 1626,
                   }

        self._testJointcalTask(2, relative_error, self.absolute_error, pa1, metrics=metrics)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
