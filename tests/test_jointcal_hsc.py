# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import inspect
import unittest
import os

from astropy import units as u

import lsst.afw.coord
import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask

import jointcalTestBase


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestHSC(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(cls.data_dir, 'hsc_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        # This value was empirically determined from the first run of jointcal on
        # this data, and will likely vary from survey to survey.
        self.dist_rms_absolute = 53e-3*u.arcsecond

        do_plot = False

        # center of the hsc validation_data catalog
        center = lsst.afw.coord.IcrsCoord(320.367492*lsst.afw.geom.degrees, 0.3131554*lsst.afw.geom.degrees)
        radius = 5*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'hsc')
        all_visits = [903334, 903336, 903338, 903342, 903344, 903346, 903986, 903988, 903990, 904010, 904014]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

    def test_jointcalTask_2_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.024
        metrics = {'collectedAstrometryRefStars': 2187,
                   'collectedPhotometryRefStars': 2187,
                   'selectedAstrometryRefStars': 2187,
                   'selectedPhotometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1151,
                   'associatedPhotometryFittedStars': 1151,
                   'selectedAstrometryFittedStars': 770,
                   'selectedPhotometryFittedStars': 770,
                   'selectedAstrometryCcdImageList': 6,
                   'selectedPhotometryCcdImageList': 6,
                   'astrometryFinalChi2': 691.12,
                   'astrometryFinalNdof': 1858,
                   'photometryFinalChi2': 3753.82,
                   'photometryFinalNdof': 504
                   }
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_11_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.134
        metrics = {'collectedAstrometryRefStars': 3649,
                   'collectedPhotometryRefStars': 3649,
                   'selectedAstrometryRefStars': 3649,
                   'selectedPhotometryRefStars': 3649,
                   'associatedAstrometryFittedStars': 2908,
                   'associatedPhotometryFittedStars': 2908,
                   'selectedAstrometryFittedStars': 2203,
                   'selectedPhotometryFittedStars': 2203,
                   'selectedAstrometryCcdImageList': 33,
                   'selectedPhotometryCcdImageList': 33,
                   'astrometryFinalChi2': 7929.656,
                   'astrometryFinalNdof': 14262,
                   'photometryFinalChi2': 16773556.5,
                   'photometryFinalNdof': 6569
                   }
        self._testJointcalTask(11, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def testJointcalTask_2_visits_no_astrometry(self):
        """Test turning off fitting astrometry."""
        pa1 = 0.024
        metrics = {'collectedPhotometryRefStars': 2187,
                   'selectedPhotometryRefStars': 2187,
                   'associatedPhotometryFittedStars': 1151,
                   'selectedPhotometryFittedStars': 770,
                   'selectedPhotometryCcdImageList': 6,
                   'photometryFinalChi2': 3753.82,
                   'photometryFinalNdof': 504
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
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collectedAstrometryRefStars': 2187,
                   'selectedAstrometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1151,
                   'selectedAstrometryFittedStars': 770,
                   'selectedAstrometryCcdImageList': 6,
                   'astrometryFinalChi2': 691.1210,
                   'astrometryFinalNdof': 1858,
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

    def test_jointcalTask_2_visits_gaia_refcat(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryRefObjLoader.retarget(LoadIndexedReferenceObjectsTask)

        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-config.py')
        self.other_args.extend(['--configfile', test_config])
        dist_rms_relative = 17e-3*u.arcsecond
        # NOTE: PA1 is slightly different here, because the number of SDSS
        # cross-matches within 0.1" goes down after we apply the GAIA-fit WCS.
        pa1 = 0.02405
        metrics = {'collectedAstrometryRefStars': 1425,
                   'collectedPhotometryRefStars': 2187,
                   'selectedAstrometryRefStars': 1425,
                   'selectedPhotometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1151,
                   'associatedPhotometryFittedStars': 1151,
                   'selectedAstrometryFittedStars': 645,
                   'selectedPhotometryFittedStars': 770,
                   'selectedAstrometryCcdImageList': 6,
                   'selectedPhotometryCcdImageList': 6,
                   'astrometryFinalChi2': 435.01995,
                   'astrometryFinalNdof': 1412,
                   'photometryFinalChi2': 3753.82,
                   'photometryFinalNdof': 504
                   }
        # NOTE: The astrometry/photometry tests are computed using the a.net SDSS refcat,
        # so the absolute astrometry RMS will be larger (because GAIA is better, so
        # comparing against SDSS will look "worse").
        dist_rms_absolute = 56e-3*u.arcsecond
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_2_visits_no_photometry_match_cut_10(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.matchCut = 10.0  # TODO: once DM-6885 is fixed, we need to put `*lsst.afw.geom.arcseconds`
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collectedAstrometryRefStars': 2187,
                   'selectedAstrometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1151,
                   'selectedAstrometryFittedStars': 790,
                   'selectedAstrometryCcdImageList': 6,
                   'astrometryFinalChi2': 690.509,
                   'astrometryFinalNdof': 1856,
                   }
        pa1 = None
        self._testJointcalTask(2, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_no_photometry(self):
        """3 visit, default config to compare with min_measurements_3 test."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.minMeasurements = 2
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collectedAstrometryRefStars': 2187,
                   'selectedAstrometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1270,
                   'selectedAstrometryFittedStars': 946,
                   'selectedAstrometryCcdImageList': 8,
                   'astrometryFinalChi2': 1229.212,
                   'astrometryFinalNdof': 3008,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)

    def test_jointcalTask_3_visits_no_photometry_min_measurements_3(self):
        """Raising min_measurements to 3 will reduce the number of selected
        fitted stars (and thus the chisq and Ndof), but should not change the
        other values."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.minMeasurements = 3
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collectedAstrometryRefStars': 2187,
                   'selectedAstrometryRefStars': 2187,
                   'associatedAstrometryFittedStars': 1270,
                   'selectedAstrometryFittedStars': 696,
                   'selectedAstrometryCcdImageList': 8,
                   'astrometryFinalChi2': 1047.57,
                   'astrometryFinalNdof': 2526,
                   }
        pa1 = None
        self._testJointcalTask(3, dist_rms_relative, self.dist_rms_absolute, pa1, metrics=metrics)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
