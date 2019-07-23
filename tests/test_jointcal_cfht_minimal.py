# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Test with a minimal catalog extracted from cfht."""
import inspect
import unittest
import os

import lsst.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase
from lsst.jointcal import jointcal
import lsst.jointcal.testUtils


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


@unittest.skipUnless(lsst.jointcal.testUtils.canRunTests(), "obs_cfht not available to use cfht_minimal.")
class JointcalTestCFHTMinimal(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    """
    Test with a stripped down CFHT dataset containing 3 stars, so by-hand
    calculation of metrics is possible.

    See `notebooks/cfht_minimal_direct_calculation.ipynb` for numpy-based
    computations of chi2, etc. using this dataset.
    """
    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/data')

    def setUp(self):
        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.geom.SpherePoint(214.884832, 52.6622199, lsst.geom.degrees)
        radius = 3*lsst.geom.degrees

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
        self.config.photometryModel = "simpleFlux"
        self.config.doAstrometry = False
        self.config.writeInitMatrix = True  # write Hessian/gradient files
        self.jointcalStatistics.do_astrometry = False

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        metrics = {'collected_photometry_refStars': 346,
                   'selected_photometry_refStars': 2,
                   'associated_photometry_fittedStars': 2,
                   'selected_photometry_fittedStars': 2,
                   'selected_photometry_ccdImages': 2,
                   'photometry_final_chi2': 2.336915,
                   'photometry_final_ndof': 1
                   }

        # The output repo is named after this method.
        caller = inspect.stack()[0].function
        # we use _runJointcalTask instead of _test here because we aren't doing
        # full calulation of PA1: the above chi2 is exact.
        self._runJointcalTask(2, caller, metrics=metrics)

        # Check that the Hessian/gradient files were written.
        self.assertTrue(os.path.exists("photometry_preinit-mat.txt"))
        os.remove("photometry_preinit-mat.txt")
        self.assertTrue(os.path.exists("photometry_preinit-grad.txt"))
        os.remove("photometry_preinit-grad.txt")
        self.assertTrue(os.path.exists("photometry_postinit-mat.txt"))
        os.remove("photometry_postinit-mat.txt")
        self.assertTrue(os.path.exists("photometry_postinit-grad.txt"))
        os.remove("photometry_postinit-grad.txt")

        # Check that the config was persisted, we can read it, and it matches the settings above
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        self.assertTrue(os.path.exists(os.path.join(output_dir, 'config/jointcal.py')))
        butler = lsst.daf.persistence.Butler(output_dir)
        config = butler.get('jointcal_config')
        self.assertEqual(config.photometryModel, self.config.photometryModel)
        self.assertEqual(config.doAstrometry, self.config.doAstrometry)
        self.assertEqual(config.writeInitMatrix, self.config.writeInitMatrix)

    def test_jointcalTask_2_visits_photometry_magnitude(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.photometryModel = "simpleMagnitude"
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        metrics = {'collected_photometry_refStars': 346,
                   'selected_photometry_refStars': 2,
                   'associated_photometry_fittedStars': 2,
                   'selected_photometry_fittedStars': 2,
                   'selected_photometry_ccdImages': 2,
                   'photometry_final_chi2': 2.23008,
                   'photometry_final_ndof': 1
                   }

        # The output repo is named after this method.
        caller = inspect.stack()[0].function
        self._runJointcalTask(2, caller, metrics=metrics)

    def test_jointcalTask_fails_raise(self):
        """Raise an exception if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryModel = "simpleFlux"
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # The output repo is named after this method.
        caller = inspect.stack()[0].function
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/config.py')
        self.configfiles = [test_config] + self.configfiles
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config', '--configfile', *self.configfiles,
                '--doraise',
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        with self.assertRaises(RuntimeError):
            jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)

    def test_jointcalTask_fails_no_raise(self):
        """exitStatus=1 if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryModel = "simpleFlux"
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # The output repo is named after this method.
        caller = inspect.stack()[0].function
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/config.py')
        self.configfiles = [test_config] + self.configfiles
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config', '--configfile', *self.configfiles,
                '--noExit',  # have to specify noExit, otherwise the test quits
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)
        self.assertEqual(result.resultList[0].exitStatus, 1)

    def test_jointcalTask_fails_no_raise_no_return_results(self):
        """exitStatus=1 if there is no data to process."""
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.setDefaults()
        self.config.photometryModel = "simpleFlux"
        self.config.sourceSelector['astrometry'].minSnr = 10000
        self.config.doAstrometry = False

        # The output repo is named after this method.
        caller = inspect.stack()[0].function
        nCatalogs = 2
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/config.py')
        self.configfiles = [test_config] + self.configfiles
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config', '--configfile', *self.configfiles,
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
