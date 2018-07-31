# See COPYRIGHT file at the top of the source tree.
import unittest
import unittest.mock

import lsst.log
import lsst.utils

import lsst.jointcal
from lsst.jointcal import MinimizeResult
import lsst.jointcal.chi2
import lsst.jointcal.testUtils


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class TestJointcalIterateFit(lsst.utils.tests.TestCase):
    def setUp(self):
        struct = lsst.jointcal.testUtils.createTwoFakeCcdImages(100, 100)
        self.ccdImageList = struct.ccdImageList
        # so that countStars() returns nonzero results
        for ccdImage in self.ccdImageList:
            ccdImage.resetCatalogForFit()

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        # disable both, so it doesn't configure any refObjLoaders
        self.config.doAstrometry = False
        self.config.doPhotometry = False
        self.jointcal = lsst.jointcal.JointcalTask(config=self.config)

        self.goodChi2 = lsst.jointcal.chi2.Chi2Statistic()
        # chi2/ndof == 2.0 should be non-bad
        self.goodChi2.chi2 = 200.0
        self.goodChi2.ndof = 100

        self.badChi2 = lsst.jointcal.chi2.Chi2Statistic()
        self.badChi2.chi2 = 600.0
        self.badChi2.ndof = 100

        self.maxSteps = 20
        self.name = "testing"
        self.whatToFit = ""  # unneeded, since we're mocking the fitter

        # Mock the fitter and association manager, so we can force particular
        # return values/exceptions. Default to "good" return values.
        self.fitter = unittest.mock.Mock(spec=lsst.jointcal.PhotometryFit)
        self.fitter.computeChi2.return_value = self.goodChi2
        self.fitter.minimize.return_value = MinimizeResult.Converged
        self.associations = unittest.mock.Mock(spec=lsst.jointcal.Associations)
        self.associations.getCcdImageList.return_value = self.ccdImageList

    def test_iterateFit_success(self):
        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.goodChi2)
        # Once for the for loop, the second time for the rank update.
        self.assertEqual(self.fitter.minimize.call_count, 2)

    def test_iterateFit_failed(self):
        self.fitter.minimize.return_value = MinimizeResult.Failed

        with self.assertRaises(RuntimeError):
            self.jointcal._iterate_fit(self.associations, self.fitter,
                                       self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(self.fitter.minimize.call_count, 1)

    def test_iterateFit_badFinalChi2(self):
        log = unittest.mock.Mock(spec=lsst.log.Log)
        self.jointcal.log = log
        self.fitter.computeChi2.return_value = self.badChi2

        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          self.maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.badChi2)
        log.info.assert_called_with("Fit completed with: %s", str(self.badChi2))
        log.error.assert_called_with("Potentially bad fit: High chi-squared/ndof.")

    def test_iterateFit_exceedMaxSteps(self):
        log = unittest.mock.Mock(spec=lsst.log.Log)
        self.jointcal.log = log
        self.fitter.minimize.return_value = MinimizeResult.Chi2Increased
        maxSteps = 3

        chi2 = self.jointcal._iterate_fit(self.associations, self.fitter,
                                          maxSteps, self.name, self.whatToFit)
        self.assertEqual(chi2, self.goodChi2)
        self.assertEqual(self.fitter.minimize.call_count, maxSteps)
        log.error.assert_called_with("testing failed to converge after %s steps" % maxSteps)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
