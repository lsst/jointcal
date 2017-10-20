# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

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
        self.dist_rms_absolute = 48.6e-3*u.arcsecond

        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.coord.IcrsCoord(214.884832*lsst.afw.geom.degrees, 52.6622199*lsst.afw.geom.degrees)
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

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 25e-3*u.arcsecond
        pa1 = 0.014
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
                   'photometryFinalChi2': 2820.84,
                   'photometryFinalNdof': 1388
                   }

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPoly(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.astrometryModel = "constrainedPoly"
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 12e-3*u.arcsecond
        pa1 = None
        metrics = {'collectedAstrometryRefStars': 825,
                   'selectedAstrometryRefStars': 825,
                   'associatedAstrometryFittedStars': 2269,
                   'selectedAstrometryFittedStars': 1239,
                   'selectedAstrometryCcdImageList': 12,
                   'astrometryFinalChi2': 1241.6,
                   'astrometryFinalNdof': 2640,
                   }

        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_constrainedPhotometry_no_astrometry(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.photometryModel = "constrained"
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        pa1 = 0.017
        metrics = {'collectedPhotometryRefStars': 825,
                   'selectedPhotometryRefStars': 825,
                   'associatedPhotometryFittedStars': 2269,
                   'selectedPhotometryFittedStars': 1239,
                   'selectedPhotometryCcdImageList': 12,
                   'photometryFinalChi2': 2640.74,
                   'photometryFinalNdof': 1328
                   }

        self._testJointcalTask(2, None, None, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
