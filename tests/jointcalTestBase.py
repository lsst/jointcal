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

__all__ = ["importRepository", "JointcalTestBase"]

import os
import shutil
import tempfile
import numpy as np

import lsst.afw.image.utils
from lsst.ctrl.mpexec import SimplePipelineExecutor
import lsst.daf.butler
from lsst.daf.butler.script import ingest_files
import lsst.pipe.base
import lsst.obs.base
import lsst.geom
from lsst.verify.bin.jobReporter import JobReporter

import lsst.jointcal


def importRepository(instrument, exportPath, exportFile, outputDir=None,
                     refcats=None, refcatPath=""):
    """Import a gen3 test repository into a test path.

    Parameters
    ----------
    instrument : `str`
        Full string name for the instrument.
    exportPath : `str`
        Path to location of repository to export.
    exportFile : `str`
        Filename of YAML butler export data, describing the data to import.
    outputDir : `str`, or `None`
        Root path to put the test repo in; appended with `testrepo/`.
        If not supplied, a temporary path is generated with mkdtemp; this path
        must be deleted by the calling code, it is not automatically removed.
    refcats : `dict` [`str`, `str`]
        Mapping of refcat name to relative path to .ecsv ingest file.
    refcatPath : `str`
        Path prefix for refcat ingest files and files linked therein.

    Returns
    -------
    repopath : `str`
        The path to the newly created butler repo.
    """
    if outputDir is None:
        repopath = tempfile.mkdtemp()
    else:
        repopath = os.path.join(outputDir, "testrepo")

    # Make the repo and retrieve a writeable Butler
    _ = lsst.daf.butler.Butler.makeRepo(repopath)
    butler = lsst.daf.butler.Butler(repopath, writeable=True)

    instrInstance = lsst.pipe.base.Instrument.from_string(instrument)
    instrInstance.register(butler.registry)

    # Register refcats first, so the `refcats` collection will exist.
    if refcats is not None:
        for name, file in refcats.items():
            graph = butler.registry.dimensions.extract(["htm7"])
            datasetType = lsst.daf.butler.DatasetType(name, graph, "SimpleCatalog",
                                                      universe=butler.registry.dimensions)
            butler.registry.registerDatasetType(datasetType)
            ingest_files(repopath, name, "refcats", file, prefix=refcatPath)
        # New butler to get the refcats collections to be visible
        butler = lsst.daf.butler.Butler(repopath, writeable=True)

    butler.import_(directory=exportPath, filename=exportFile,
                   transfer='symlink',
                   skip_dimensions={'instrument', 'detector', 'physical_filter'})
    return repopath


class JointcalTestBase:
    """
    Base class for jointcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """

    def set_output_dir(self):
        """Set the output directory to the name of the test method, in .test/
        """
        self.output_dir = os.path.join('.test', self.__class__.__name__, self.id().split('.')[-1])

    def setUp_base(self,
                   instrumentClass,
                   instrumentName,
                   input_dir="",
                   all_visits=None,
                   log_level=None,
                   where="",
                   refcats=None,
                   refcatPath="",
                   inputCollections=None,
                   outputDataId=None):
        """
        Call from your child classes's setUp() to get the necessary variables built.

        Parameters
        ----------
        instrumentClass : `str`
            The full module name of the instrument to be registered in the
            new repo. For example, "lsst.obs.subaru.HyperSuprimeCam"
        instrumentName : `str`
            The name of the instrument as it appears in the repo collections.
            For example, "HSC".
        input_dir : `str`
            Directory of input butler repository.
        all_visits : `list` [`int`]
            List of the available visits to generate the dataset query from.
        log_level : `str`
            Set to the default log level you want jointcal to produce while the
            tests are running. See the developer docs about logging for valid
            levels: https://developer.lsst.io/coding/logging.html
        where : `str`
            Data ID query for pipetask specifying the data to run on.
        refcats : `dict` [`str`, `str`]
            Mapping of refcat name to relative path to .ecsv ingest file.
        refcatPath : `str`
            Path prefix for refcat ingest files and files linked therein.
        inputCollections : `list` [`str`]
            String to use for "-i" input collections (comma delimited).
            For example, "refcats,HSC/runs/tests,HSC/calib"
        outputDataId : `dict`
            Partial dataIds for testing whether the output files were written.
        """
        self.path = os.path.dirname(__file__)

        self.instrumentClass = instrumentClass
        self.instrumentName = instrumentName

        self.input_dir = input_dir
        self.all_visits = all_visits
        self.log_level = log_level

        self.configfiles = []

        # Make unittest output more verbose failure messages to assert failures.
        self.longMessage = True

        self.where = where
        self.refcats = refcats
        self.refcatPath = refcatPath
        self.inputCollections = inputCollections

        self.outputDataId = outputDataId

        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination from other instruments.
        lsst.obs.base.FilterDefinitionCollection.reset()

        self.set_output_dir()

    def tearDown(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)

    def _test_metrics(self, result, expect):
        """Test a dictionary of "metrics" against those returned by jointcal.py

        Parameters
        ----------
        result : `dict`
            Result metric dictionary from jointcal.py
        expect : `dict`
            Expected metric dictionary; set a value to None to not test it.
        """
        for key in result:
            with self.subTest(key.metric):
                if expect[key.metric] is not None:
                    value = result[key].quantity.value
                    if isinstance(value, float):
                        self.assertFloatsAlmostEqual(value, expect[key.metric], msg=key.metric, rtol=1e-5)
                    else:
                        self.assertEqual(value, expect[key.metric], msg=key.metric)

    def _runJointcal(self, repo,
                     inputCollections, outputCollection,
                     configFiles=None, configOptions=None,
                     registerDatasetTypes=False, whereSuffix=None,
                     nJobs=1):
        """Run a pipeline via the SimplePipelineExecutor.

        Parameters
        ----------
        repo : `str`
            Gen3 Butler repository to read from/write to.
        inputCollections : `list` [`str`]
            String to use for "-i" input collections (comma delimited).
            For example, "refcats,HSC/runs/tests,HSC/calib"
        outputCollection : `str`
            Name of the output collection to write to. For example,
            "HSC/testdata/jointcal"
        configFiles : `list` [`str`], optional
            List of jointcal config files to use.
        configOptions : `dict` [`str`], optional
            Individual jointcal config options (field: value) to override.
        registerDatasetTypes : bool, optional
            Set "--register-dataset-types" when running the pipeline.
        whereSuffix : `str`, optional
            Additional parameters to the ``where`` pipetask statement.
        nJobs : `int`, optional
            Number of quanta expected to be run.

        Returns
        -------
        job : `lsst.verify.Job`
            Job containing the metric measurements from this test run.
        """
        config = lsst.jointcal.JointcalConfig()
        if configFiles:
            for file in configFiles:
                config.load(file)
        if configOptions:
            for key, value in configOptions.items():
                setattr(config, key, value)
        where = ' '.join((self.where, whereSuffix)) if whereSuffix is not None else self.where
        lsst.daf.butler.cli.cliLog.CliLog.initLog(False)
        butler = SimplePipelineExecutor.prep_butler(repo,
                                                    inputs=inputCollections,
                                                    output=outputCollection,
                                                    output_run=outputCollection.replace("all", "jointcal"))
        executor = SimplePipelineExecutor.from_task_class(lsst.jointcal.JointcalTask,
                                                          config=config,
                                                          where=where,
                                                          butler=butler)
        executor.run(register_dataset_types=registerDatasetTypes)
        # JobReporter bundles all metrics in the collection into one job.
        jobs = JobReporter(repo, outputCollection, "jointcal", "", "jointcal").run()
        # should only ever get one job output in tests, unless specified
        self.assertEqual(len(jobs), nJobs)
        # Sort the jobs, as QuantumGraph ordering is not guaranteed, and some
        # tests check metrics values for one job while producing two.
        sorted_jobs = {key: jobs[key] for key in sorted(list(jobs.keys()))}
        return list(sorted_jobs.values())[0]

    def _runJointcalTest(self, astrometryOutputs=None, photometryOutputs=None,
                         configFiles=None, configOptions=None, whereSuffix=None,
                         metrics=None, nJobs=1):
        """Create a Butler repo and run jointcal on it.

        Parameters
        ----------
        atsrometryOutputs, photometryOutputs : `dict` [`int`, `list`]
            visit: [detectors] dictionary to test that output files were saved.
            Default None means that there should be no output for that
            respective dataset type.
        configFiles : `list` [`str`], optional
            List of jointcal config files to use.
        configOptions : `dict` [`str`], optional
            Individual jointcal config options (field: value) to override.
        whereSuffix : `str`, optional
            Additional parameters to the ``where`` pipetask statement.
        metrics : `dict`, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.
        nJobs : `int`, optional
            Number of quanta expected to be run.

        Returns
        -------
        repopath : `str`
            The path to the newly created butler repo.
        """
        repopath = importRepository(self.instrumentClass,
                                    os.path.join(self.input_dir, "repo"),
                                    os.path.join(self.input_dir, "exports.yaml"),
                                    self.output_dir,
                                    refcats=self.refcats,
                                    refcatPath=self.refcatPath)
        configs = [os.path.join(self.path, "config/config.py")]
        configs.extend(self.configfiles or [])
        configs.extend(configFiles or [])
        collection = f"{self.instrumentName}/tests/all"
        job = self._runJointcal(repopath,
                                self.inputCollections,
                                collection,
                                configFiles=configs,
                                configOptions=configOptions,
                                registerDatasetTypes=True,
                                whereSuffix=whereSuffix,
                                nJobs=nJobs)

        if metrics:
            self._test_metrics(job.measurements, metrics)

        butler = lsst.daf.butler.Butler(repopath, collections=collection)
        if astrometryOutputs is not None:
            for visit, detectors in astrometryOutputs.items():
                self._check_astrometry_output(butler, visit, detectors)
        else:
            self.assertEqual(len(list(butler.registry.queryDatasets("jointcalSkyWcsCatalog"))), 0)
        if photometryOutputs is not None:
            for visit, detectors in photometryOutputs.items():
                self._check_photometry_output(butler, visit, detectors)
        else:
            self.assertEqual(len(list(butler.registry.queryDatasets("jointcalPhotoCalibCatalog"))), 0)

        return repopath

    def _check_astrometry_output(self, butler, visit, detectors):
        """Check that the WCSs were written for each detector, and only for the
        correct detectors.
        """
        dataId = self.outputDataId.copy()
        dataId['visit'] = visit

        catalog = butler.get('jointcalSkyWcsCatalog', dataId)
        for record in catalog:
            self.assertIsInstance(record.getWcs(), lsst.afw.geom.SkyWcs,
                                  msg=f"visit {visit}: {record}")
        np.testing.assert_array_equal(catalog['id'], detectors)

    def _check_photometry_output(self, butler, visit, detectors):
        """Check that the PhotoCalibs were written for each detector, and only
        for the correct detectors.
        """
        dataId = self.outputDataId.copy()
        dataId['visit'] = visit

        catalog = butler.get('jointcalPhotoCalibCatalog', dataId)
        for record in catalog:
            self.assertIsInstance(record.getPhotoCalib(), lsst.afw.image.PhotoCalib,
                                  msg=f"visit {visit}: {record}")
        np.testing.assert_array_equal(catalog['id'], detectors)
