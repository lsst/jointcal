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
import unittest
import os

import lsst.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase
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

        where = "instrument='MegaPrime' and tract=0 and skymap='discrete'"  # and detector=12"
        inputCollections = ["singleFrame", "skymaps"]
        refcats = {"sdss_dr9_fink_v5b": os.path.join(input_dir, "sdss-dr9-fink-v5b.ecsv")}

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot,
                        inputCollections=inputCollections,
                        refcats=refcats,
                        where=where,
                        refcatPath=input_dir,
                        log_level="debug")
        test_config = os.path.join(os.path.dirname(__file__), 'config/cfht_minimal-config.py')
        self.configfiles.append(test_config)

    def test_jointcalTask_2_visits_photometry(self):
        configOptions = {"doAstrometry": False, "photometryModel": "simpleFlux",
                         "writeInitMatrix": True}

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        # This chi2 is exact, from calculations in cfht_minimal_direct_calculation.ipynb
        metrics = {'photometry_collected_refStars': 346,
                   'photometry_prepared_refStars': 2,
                   'photometry_matched_fittedStars': 2,
                   'photometry_prepared_fittedStars': 2,
                   'photometry_prepared_ccdImages': 2,
                   'photometry_final_chi2': 2.336915,
                   'photometry_final_ndof': 1
                   }

        repo = self._runGen3Jointcal("lsst.obs.cfht.MegaPrime", "MegaPrime",
                                     configOptions=configOptions, metrics=metrics)

        # Check that the Hessian/gradient files were written.
        self.assertTrue(os.path.exists("photometry_preinit-0_r.MP9601-mat.txt"))
        os.remove("photometry_preinit-0_r.MP9601-mat.txt")
        self.assertTrue(os.path.exists("photometry_preinit-0_r.MP9601-grad.txt"))
        os.remove("photometry_preinit-0_r.MP9601-grad.txt")
        self.assertTrue(os.path.exists("photometry_postinit-0_r.MP9601-mat.txt"))
        os.remove("photometry_postinit-0_r.MP9601-mat.txt")
        self.assertTrue(os.path.exists("photometry_postinit-0_r.MP9601-grad.txt"))
        os.remove("photometry_postinit-0_r.MP9601-grad.txt")

        # Check that the config was persisted, and that it matches the settings above
        config = lsst.jointcal.jointcal.JointcalConfig()
        for key, value in configOptions.items():
            setattr(config, key, value)
        butler = lsst.daf.butler.Butler(repo, collections=['MegaPrime/tests/all'])
        output_config = butler.get('jointcal_config')
        self.assertEqual(output_config.photometryModel, config.photometryModel)
        self.assertEqual(output_config.doAstrometry, config.doAstrometry)
        self.assertEqual(output_config.writeInitMatrix, config.writeInitMatrix)

    def test_jointcalTask_2_visits_photometry_magnitude(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        configOptions = {"doAstrometry": False, "photometryModel": "simpleMagnitude"}

        # NOTE: ndof==1 from 4 fit parameters (2 model, 2 fittedStar), and
        # 5 degrees-of-freedom (3 star measurements, with 2 reference stars).
        metrics = {'photometry_collected_refStars': 346,
                   'photometry_prepared_refStars': 2,
                   'photometry_matched_fittedStars': 2,
                   'photometry_prepared_fittedStars': 2,
                   'photometry_prepared_ccdImages': 2,
                   'photometry_final_chi2': 2.23008,
                   'photometry_final_ndof': 1
                   }

        self._runGen3Jointcal("lsst.obs.cfht.MegaPrime", "MegaPrime",
                              configOptions=configOptions, metrics=metrics)

    def test_jointcalTask_fails_raise(self):
        """Raise an exception if there is no data to process."""
        self.configfiles.append(os.path.join(self.path, 'config/minSnr.py'))
        with self.assertRaisesRegex(RuntimeError, "No data to process"):
            self._runGen3Jointcal("lsst.obs.cfht.MegaPrime", "MegaPrime")


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
