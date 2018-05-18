"""Test with a minimal catalog extracted from cfht."""
import inspect
import unittest
import os
import contextlib

import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask

import jointcalTestBase
from lsst.jointcal import jointcal


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHTMinimal(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/data')
            anet_data_dir = lsst.utils.getPackageDir('testdata_jointcal')
            os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(anet_data_dir, 'cfht_and_index')
        except lsst.pex.exceptions.NotFoundError:
            raise unittest.SkipTest("testdata_jointcal not setup")

    def setUp(self):
        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.geom.SpherePoint(214.884832, 52.6622199, lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(self.data_dir, 'cfht_minimal')
        all_visits = [849375, 850587]
        other_args = ['ccd=12']

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        other_args=other_args,
                        do_plot=do_plot, log_level="debug")

    def tearDown(self):
        super()
        with contextlib.suppress(FileNotFoundError):
            os.remove("photometry_preinit-mat.txt")
            os.remove("photometry_preinit-grad.txt")
            os.remove("photometry_postinit-mat.txt")
            os.remove("photometry_postinit-grad.txt")

    def test_jointcalTask_2_visits_photometry(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.doAstrometry = False
        self.config.writeInitMatrix = True  # write Hessian/gradient files
        self.jointcalStatistics.do_astrometry = False

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        metrics = {'collected_photometry_refStars': 183,
                   'selected_photometry_refStars': 2,
                   'associated_photometry_fittedStars': 2,
                   'selected_photometry_fittedStars': 2,
                   'selected_photometry_ccdImages': 2,
                   'photometry_final_chi2': 2.336915,
                   'photometry_final_ndof': 1
                   }

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        self._runJointcalTask(2, caller, metrics=metrics)

        # Check that the Hessian/gradient files were written.
        self.assertTrue(os.path.exists("photometry_preinit-mat.txt"))
        self.assertTrue(os.path.exists("photometry_preinit-grad.txt"))
        self.assertTrue(os.path.exists("photometry_postinit-mat.txt"))
        self.assertTrue(os.path.exists("photometry_postinit-grad.txt"))

    def test_jointcalTask_fails_raise(self):
        """Raise an exception if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config',
                '--doraise',
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        with self.assertRaises(RuntimeError):
            jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)

    def test_jointcalTask_fails_no_raise(self):
        """exitStatus=1 if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config',
                '--noExit',  # have to specify noExit, otherwise the test quits
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)
        self.assertEqual(result.resultList[0].exitStatus, 1)

    def test_jointcalTask_fails_no_raise_no_return_results(self):
        """exitStatus=1 if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryRefObjLoader.retarget(LoadAstrometryNetObjectsTask)
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config',
                '--noExit',  # have to specify noExit, otherwise the test quits
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, config=self.config)
        self.assertEqual(result.resultList[0].exitStatus, 1)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
