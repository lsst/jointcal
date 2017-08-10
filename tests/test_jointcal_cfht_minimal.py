"""Test with a minimal catalog extracted from cfht."""
from __future__ import division, absolute_import, print_function

import inspect
import unittest
import os

import lsst.afw.coord
import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHTMinimal(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'cfht_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.coord.IcrsCoord(214.884832*lsst.afw.geom.degrees, 52.6622199*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'cfht_minimal')
        all_visits = [849375, 850587]
        other_args = ['ccd=12']

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot, log_level="debug")

    def test_jointcalTask_2_visits_photometry(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        metrics = {'collectedPhotometryRefStars': 183,
                   'selectedPhotometryRefStars': 183,
                   'associatedPhotometryFittedStars': 2,
                   'selectedPhotometryFittedStars': 2,
                   'selectedPhotometryCcdImageList': 2,
                   'photometryFinalChi2': 2.336915,
                   'photometryFinalNdof': 1
                   }

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        self._runJointcalTask(2, caller, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
