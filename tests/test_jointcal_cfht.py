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
        # This value was empirically determined from the first run of jointcal on
        # this data, and will likely vary from survey to survey.
        self.dist_rms_absolute = 48e-3*u.arcsecond

        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.coord.IcrsCoord(214.884832*lsst.afw.geom.degrees, 52.6622199*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'cfht')
        all_visits = [849375, 850587]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

    def test_jointcalTask_2_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 25e-3*u.arcsecond
        pa1 = 0.019
        metrics = {'collectedAstrometryRefStars': 825,
                   'collectedPhotometryRefStars': 825,
                   'selectedAstrometryRefStars': 825,
                   'selectedPhotometryRefStars': 825,
                   'associatedAstrometryFittedStars': 2269,
                   'associatedPhotometryFittedStars': 2269,
                   'selectedAstrometryFittedStars': 1239,
                   'selectedPhotometryFittedStars': 1239,
                   'selectedAstrometryCcdImageList': 12,
                   'selectedPhotometryCcdImageList': 12,
                   'astrometryFinalChi2': 1150.62,
                   'astrometryFinalNdof': 2550,
                   'photometryFinalChi2': 13363.6,
                   'photometryFinalNdof': 1089
                   }

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
