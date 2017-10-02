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


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestLSSTSim(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'twinkles1_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # We don't want the absolute astrometry to become significantly worse
        # than the single-epoch astrometry (about 0.040").
        # This value was empirically determined from the first run of jointcal on
        # this data, and will likely vary from survey to survey.
        self.dist_rms_absolute = 42e-3*u.arcsecond

        do_plot = False

        # position of the Twinkles run 1 catalog
        center = lsst.afw.coord.IcrsCoord(53.00914*lsst.afw.geom.degrees, -27.43895*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'twinkles1')
        all_visits = list(range(840, 850))
        other_args = ['raft=2,2', 'sensor=1,1', 'filter=r']

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot)

    @unittest.skip('jointcal currently fails (may segfault) if only given one catalog!')
    def testJointcalTask_1_visits(self):
        dist_rms_relative = 0*u.arcsecond  # there is no such thing as a "relative" test for 1 catalog.
        pa1 = 2.64e-3
        self._testJointcalTask(1, dist_rms_relative, self.dist_rms_absolute, pa1)

    def testJointcalTask_2_visits(self):
        # NOTE: The relative RMS limits were empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 9.7e-3*u.arcsecond
        # TODO: reinstate photometry tests once DM-11397 is fixed.
        # NOTE: the measured values of the metrics may change with the new fitter.
        pa1 = None
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False
        # pa1 = 2.64e-3
        # 'collectedPhotometryRefStars': 1686,
        # 'selectedPhotometryRefStars': 1686,
        # 'associatedPhotometryFittedStars': 1365,
        # 'selectedPhotometryFittedStars': 789,
        # 'selectedPhotometryCcdImageList': 2,
        # 'photometryFinalChi2': 2310.628,
        # 'photometryFinalNdof': 727
        metrics = {'collectedAstrometryRefStars': 1686,
                   'selectedAstrometryRefStars': 1686,
                   'associatedAstrometryFittedStars': 1365,
                   'selectedAstrometryFittedStars': 789,
                   'selectedAstrometryCcdImageList': 2,
                   'astrometryFinalChi2': 742.4992,
                   'astrometryFinalNdof': 1698,
                   }
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def testJointcalTask_10_visits(self):
        dist_rms_relative = 7.9e-3*u.arcsecond
        # TODO: reinstate photometry tests once DM-11397 is fixed.
        # NOTE: the measured values of the metrics may change with the new fitter.
        pa1 = None
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False
        # pa1 = 2.64e-3
        # 'collectedPhotometryRefStars': 1686,
        # 'selectedPhotometryRefStars': 1686,
        # 'associatedPhotometryFittedStars': 1823,
        # 'selectedPhotometryFittedStars': 1506,
        # 'selectedPhotometryCcdImageList': 10,
        # 'photometryFinalChi2': 35321.947,
        # 'photometryFinalNdof': 9140
        metrics = {'collectedAstrometryRefStars': 1686,
                   'selectedAstrometryRefStars': 1686,
                   'associatedAstrometryFittedStars': 1823,
                   'selectedAstrometryFittedStars': 1506,
                   'selectedAstrometryCcdImageList': 10,
                   'astrometryFinalChi2': 7262.2075,
                   'astrometryFinalNdof': 18260,
                   }
        self._testJointcalTask(10, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    @unittest.skip("A.net reference catalog missing flux errors. Unskip once DM-.")
    def testJointcalTask_2_visits_no_astrometry(self):
        """Test turning off fitting astrometry."""
        pa1 = 2.64e-3
        metrics = {'collectedPhotometryRefStars': 1686,
                   'selectedPhotometryRefStars': 1686,
                   'associatedPhotometryFittedStars': 1365,
                   'selectedPhotometryFittedStars': 789,
                   'selectedPhotometryCcdImageList': 2,
                   'photometryFinalChi2': 2310.6280,
                   'photometryFinalNdof': 727
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller, metrics=metrics)

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

    def testJointcalTask_2_visits_no_photometry(self):
        """Test turning off fitting photometry."""
        dist_rms_relative = 9.7e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 1686,
                   'selected_astrometry_refStars': 1686,
                   'associated_astrometry_fittedStars': 1365,
                   'selected_astrometry_fittedStars': 789,
                   'selected_astrometry_ccdImages': 2,
                   'astrometry_chisq': 742.4992,
                   'astrometry_ndof': 1698
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller, metrics=metrics)
        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)

        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        self.assertLess(rms_result.dist_relative, dist_rms_relative)
        self.assertLess(rms_result.dist_absolute, self.dist_rms_absolute)
        self.assertIsNone(rms_result.pa1)

        for data_ref in data_refs:
            calib = data_ref.get('wcs').getCalib()
            blank_calib = lsst.afw.image.Calib()
            self.assertEqual(calib, blank_calib)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
