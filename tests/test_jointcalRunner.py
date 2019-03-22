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

import os.path
import tempfile
import unittest
import unittest.mock

import lsst.log
import lsst.utils

import lsst.daf.persistence
import lsst.jointcal
from lsst.meas.algorithms import DatasetConfig
import lsst.pipe.base

current_dir = os.path.abspath(os.path.dirname(__file__))


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class TestJointcalRunner(lsst.utils.tests.TestCase):
    """Test that JointcalRunner calls JointcalTask with appropriate arguments.

    This test mocks a butler to force JointcalRunner to be called with
    specific sets of parameters (e.g. 1 or 2 tracts), and mocks JointcalTask to
    check that it is called the correct number of times (i.e. once per tract).
    """
    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.join(current_dir, 'data')

    def setUp(self):
        def mock_get(arg):
            """Return things that sub-tasks need from butler.get()."""
            if arg == 'verify_job_filename':
                return 'some/path/filename.job'
            elif arg == 'ref_cat_config':
                result = unittest.mock.Mock(DatasetConfig)
                result.indexer = DatasetConfig().indexer
                return result
            else:
                return None

        def mock_dataRef(datasetType, level=None, dataId=None):
            """Return a fake dataRef that always exists."""
            ref = lsst.daf.persistence.ButlerDataRef(None, dataId)
            ref.datasetExists = lambda: True
            return ref

        def mock_subset(datasetType, dataId=None):
            """Return a fake butler subset containing one dataRef."""
            return [mock_dataRef(datasetType, dataId=dataId)]

        butlerPatcher = unittest.mock.patch("lsst.daf.persistence.Butler", autospec=True)
        self.butler = butlerPatcher.start()
        self.addCleanup(butlerPatcher.stop)
        self.butler.return_value.get = mock_get
        self.butler.return_value.subset = mock_subset
        self.butler.return_value.dataRef = mock_dataRef
        self.butler.return_value.getKeys.return_value = {'visit': int, 'ccd': int}
        # static methods: mock on the class, not the return_value
        self.butler.getMapperClass.return_value.getPackageName.return_value = 'jointcal'
        self.butler.getMapperClass.return_value.getCameraName.return_value = 'cfht'

    def prep_jointcal(self, tracts=None):
        """Prepare a jointcal mock to be called by JointcalRunner.

        We use the `tests/data/cfht_minimal` repo to provide a "real" refcat,
        and use dataIds associated with it, even though we mock the butler.

        Parameters
        ----------
        tracts : `list` [`int`]
            List of tracts to build DataRefs for.
        """
        configfile = os.path.join(current_dir, 'config/config.py')

        input_dir = os.path.join(self.data_dir, 'cfht_minimal')
        output_dir = tempfile.mkdtemp()  # we don't care about any outputs

        visits = '849375^850587'
        args = [input_dir, '--output', output_dir, '--doraise', '--configfile=%s'%configfile,
                '--id', 'visit=%s'%visits, 'ccd=12']
        if tracts is not None:
            args.append('tract=%s'%'^'.join(str(tract) for tract in tracts))

        config = lsst.jointcal.JointcalConfig()
        parser = lsst.jointcal.JointcalTask._makeArgumentParser()
        parsedArgs = parser.parse_args(config, args=args)

        jointcalPatcher = unittest.mock.patch("lsst.jointcal.JointcalTask", autospec=True)
        JointcalTask = jointcalPatcher.start()
        self.addCleanup(jointcalPatcher.stop)
        self.return_value = 100
        result = lsst.pipe.base.Struct(exitStatus=self.return_value, job=unittest.mock.Mock())
        JointcalTask.return_value.runDataRef.return_value = result
        JointcalTask.return_value.log = unittest.mock.Mock()  # jointcal.log is created in __init__

        return parsedArgs, JointcalTask

    def testJointcalRunnerOneTract(self):
        parsedArgs, JointcalTask = self.prep_jointcal(tracts=[5])
        runner = lsst.jointcal.JointcalRunner(JointcalTask, parsedArgs)
        result = runner.run(parsedArgs)
        self.assertEqual(result[0].exitStatus, self.return_value)
        JointcalTask.return_value.runDataRef.assert_called_once()

    def testJointcalRunnerTwoTracts(self):
        parsedArgs, JointcalTask = self.prep_jointcal(tracts=[5, 6])
        runner = lsst.jointcal.JointcalRunner(JointcalTask, parsedArgs)
        result = runner.run(parsedArgs)
        self.assertEqual(result[0].exitStatus, self.return_value)
        self.assertEqual(JointcalTask.return_value.runDataRef.call_count, 2)

    def testJointcalRunnerThreeTracts(self):
        parsedArgs, JointcalTask = self.prep_jointcal(tracts=[5, 6, 7])
        runner = lsst.jointcal.JointcalRunner(JointcalTask, parsedArgs)
        result = runner.run(parsedArgs)
        self.assertEqual(result[0].exitStatus, self.return_value)
        self.assertEqual(JointcalTask.return_value.runDataRef.call_count, 3)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
