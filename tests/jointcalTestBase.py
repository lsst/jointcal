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

import copy
import os
import shutil
import tempfile

import lsst.afw.image.utils
from lsst.ctrl.mpexec import SimplePipelineExecutor
import lsst.daf.butler
from lsst.daf.butler.script import ingest_files
import lsst.obs.base
import lsst.geom
from lsst.verify.bin.jobReporter import JobReporter

from lsst.jointcal import jointcal, utils


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
        repopath = os.path.join(outputDir, 'testrepo')

    # Make the repo and retrieve a writeable Butler
    _ = lsst.daf.butler.Butler.makeRepo(repopath)
    butler = lsst.daf.butler.Butler(repopath, writeable=True)

    instrInstance = lsst.obs.base.utils.getInstrument(instrument)
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
        self.output_dir = os.path.join('.test', self.__class__.__name__, self.id().split('.')[-1])

    def setUp_base(self, center, radius,
                   match_radius=0.1*lsst.geom.arcseconds,
                   input_dir="",
                   all_visits=None,
                   other_args=None,
                   do_plot=False,
                   log_level=None,
                   where="",
                   refcats=None,
                   refcatPath="",
                   inputCollections=None):
        """
        Call from your child classes's setUp() to get the necessary variables built.

        Parameters
        ----------
        center : `lsst.geom.SpherePoint`
            Center of the reference catalog.
        radius : `lsst.geom.Angle`
            Radius from center to load reference catalog objects inside.
        match_radius : `lsst.geom.Angle`
            matching radius when calculating RMS of result.
        input_dir : `str`
            Directory of input butler repository.
        all_visits : `list` [`int`]
            List of the available visits to generate the parseAndRun arguments.
        other_args : `list` [`str`]
            Optional other arguments for the butler dataId.
        do_plot : `bool`
            Set to True for a comparison plot and some diagnostic numbers.
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
        """
        self.path = os.path.dirname(__file__)

        self.center = center
        self.radius = radius
        self.jointcalStatistics = utils.JointcalStatistics(match_radius, verbose=True)
        self.input_dir = input_dir
        self.all_visits = all_visits
        if other_args is None:
            other_args = []
        self.other_args = other_args
        self.do_plot = do_plot
        self.log_level = log_level
        # Signal/Noise (flux/fluxErr) for sources to be included in the RMS cross-match.
        # 100 is a balance between good centroids and enough sources.
        self.flux_limit = 100

        # Individual tests may want to tweak the config that is passed to parseAndRun().
        self.config = None
        self.configfiles = []

        # Append `msg` arguments to assert failures.
        self.longMessage = True

        self.where = where
        self.refcats = refcats
        self.refcatPath = refcatPath
        self.inputCollections = inputCollections

        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination from other instruments.
        lsst.obs.base.FilterDefinitionCollection.reset()

        self.set_output_dir()

    def tearDown(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)

        if getattr(self, 'reference', None) is not None:
            del self.reference
        if getattr(self, 'oldWcsList', None) is not None:
            del self.oldWcsList
        if getattr(self, 'jointcalTask', None) is not None:
            del self.jointcalTask
        if getattr(self, 'jointcalStatistics', None) is not None:
            del self.jointcalStatistics
        if getattr(self, 'config', None) is not None:
            del self.config

    def _testJointcalTask(self, nCatalogs, dist_rms_relative, dist_rms_absolute, pa1,
                          metrics=None):
        """
        Test parseAndRun for jointcal on nCatalogs.

        Checks relative and absolute astrometric error (arcsec) and photometric
        repeatability (PA1 from the SRD).

        Parameters
        ----------
        nCatalogs : `int`
            Number of catalogs to run jointcal on. Used to construct the "id"
            field for parseAndRun.
        dist_rms_relative : `astropy.Quantity`
            Minimum relative astrometric rms post-jointcal to pass the test.
        dist_rms_absolute : `astropy.Quantity`
            Minimum absolute astrometric rms post-jointcal to pass the test.
        pa1 : `float`
            Minimum PA1 (from Table 14 of the Science Requirements Document:
            https://ls.st/LPM-17) post-jointcal to pass the test.
        metrics : `dict`, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.

        Returns
        -------
        dataRefs : `list` [`lsst.daf.persistence.ButlerDataRef`]
            The dataRefs that were processed.
        """

        resultFull = self._runJointcalTask(nCatalogs, metrics=metrics)
        result = resultFull.resultList[0].result  # shorten this very long thing

        def compute_statistics(refObjLoader):
            refCat = refObjLoader.loadSkyCircle(self.center,
                                                self.radius,
                                                result.defaultFilter.bandLabel,
                                                epoch=result.epoch).refCat
            rms_result = self.jointcalStatistics.compute_rms(result.dataRefs, refCat)
            # Make plots before testing, if requested, so we still get plots if tests fail.
            if self.do_plot:
                self._plotJointcalTask(result.dataRefs, result.oldWcsList)
            return rms_result

        # we now have different astrometry/photometry refcats, so have to
        # do these calculations separately
        if self.jointcalStatistics.do_astrometry:
            refObjLoader = result.astrometryRefObjLoader
            # preserve do_photometry for the next `if`
            temp = copy.copy(self.jointcalStatistics.do_photometry)
            self.jointcalStatistics.do_photometry = False
            rms_result = compute_statistics(refObjLoader)
            self.jointcalStatistics.do_photometry = temp  # restore do_photometry

            if dist_rms_relative is not None and dist_rms_absolute is not None:
                self.assertLess(rms_result.dist_relative, dist_rms_relative)
                self.assertLess(rms_result.dist_absolute, dist_rms_absolute)

        if self.jointcalStatistics.do_photometry:
            refObjLoader = result.photometryRefObjLoader
            self.jointcalStatistics.do_astrometry = False
            rms_result = compute_statistics(refObjLoader)

            if pa1 is not None:
                self.assertLess(rms_result.pa1, pa1)

        return result.dataRefs

    def _runJointcalTask(self, nCatalogs, metrics=None):
        """
        Run jointcalTask on nCatalogs, with the most basic tests.
        Tests for non-empty result list, and that the basic metrics are correct.

        Parameters
        ----------
        nCatalogs : `int`
            Number of catalogs to test on.
        metrics : `dict`, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.

        Returns
        -------
        result : `pipe.base.Struct`
            The structure returned by jointcalTask.run()
        """
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        if self.log_level is not None:
            self.other_args.extend(['--loglevel', 'jointcal=%s'%self.log_level])

        #  Place default configfile first so that specific subclass configfiles are applied after
        test_config = os.path.join(self.path, 'config/config.py')
        configfiles = [test_config] + self.configfiles

        args = [self.input_dir, '--output', self.output_dir,
                '--clobber-versions', '--clobber-config',
                '--doraise', '--configfile', *configfiles,
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)
        job = result.resultList[0].result.job
        self._test_metrics(job.measurements, metrics)

        return result

    def _plotJointcalTask(self, data_refs, oldWcsList):
        """
        Plot the results of a jointcal run.

        Parameters
        ----------
        data_refs : `list` [`lsst.daf.persistence.ButlerDataRef`]
            The dataRefs that were processed.
        oldWcsList : `list` [`lsst.afw.image.SkyWcs`]
            The original WCS from each dataRef.
        """
        plot_dir = os.path.join('.test', self.__class__.__name__, 'plots')
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)
        self.jointcalStatistics.make_plots(data_refs, oldWcsList, name=self.id(), outdir=plot_dir)
        print("Plots saved to: {}".format(plot_dir))

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
            if expect[key.metric] is not None:
                value = result[key].quantity.value
                if isinstance(value, float):
                    self.assertFloatsAlmostEqual(value, expect[key.metric], msg=key.metric, rtol=1e-5)
                else:
                    self.assertEqual(value, expect[key.metric], msg=key.metric)

    def _runPipeline(self, repo,
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
        for file in configFiles:
            config.load(file)
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
        return list(jobs.values())[0]

    def _runGen3Jointcal(self, instrumentClass, instrumentName,
                         configFiles=None, configOptions=None, whereSuffix=None,
                         metrics=None, nJobs=1):
        """Create a Butler repo and run jointcal on it.

        Parameters
        ----------
        instrumentClass : `str`
            The full module name of the instrument to be registered in the
            new repo. For example, "lsst.obs.subaru.HyperSuprimeCam"
        instrumentName : `str`
            The name of the instrument as it appears in the repo collections.
            For example, "HSC".
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
        repopath = importRepository(instrumentClass,
                                    os.path.join(self.input_dir, "repo"),
                                    os.path.join(self.input_dir, "exports.yaml"),
                                    self.output_dir,
                                    refcats=self.refcats,
                                    refcatPath=self.refcatPath)
        configs = [os.path.join(self.path, "config/config-gen3.py")]
        configs.extend(self.configfiles or [])
        configs.extend(configFiles or [])
        job = self._runPipeline(repopath,
                                self.inputCollections,
                                f"{instrumentName}/tests/all",
                                configFiles=configs,
                                configOptions=configOptions,
                                registerDatasetTypes=True,
                                whereSuffix=whereSuffix,
                                nJobs=nJobs)

        if metrics:
            self._test_metrics(job.measurements, metrics)

        return repopath
