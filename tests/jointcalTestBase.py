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

import copy
import os
import shutil

import click.testing

import lsst.afw.image.utils
import lsst.ctrl.mpexec.cli.pipetask
import lsst.daf.butler
import lsst.obs.base
import lsst.geom

from lsst.jointcal import jointcal, utils


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
                   log_level=None):
        """
        Call from your child classes's setUp() to get the necessary variables built.

        Parameters
        ----------
        center : lsst.geom.SpherePoint
            Center of the reference catalog.
        radius : lsst.geom.Angle
            Radius from center to load reference catalog objects inside.
        match_radius : lsst.geom.Angle
            matching radius when calculating RMS of result.
        input_dir : str
            Directory of input butler repository.
        all_visits : list
            List of the available visits to generate the parseAndRun arguments.
        other_args : list
            Optional other arguments for the butler dataId.
        do_plot : bool
            Set to True for a comparison plot and some diagnostic numbers.
        log_level : str
            Set to the default log level you want jointcal to produce while the
            tests are running. See the developer docs about logging for valid
            levels: https://developer.lsst.io/coding/logging.html
        """
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
        nCatalogs : int
            Number of catalogs to run jointcal on. Used to construct the "id"
            field for parseAndRun.
        dist_rms_relative : astropy.Quantity
            Minimum relative astrometric rms post-jointcal to pass the test.
        dist_rms_absolute : astropy.Quantity
            Minimum absolute astrometric rms post-jointcal to pass the test.
        pa1 : float
            Minimum PA1 (from Table 14 of the Science Requirements Document:
            https://ls.st/LPM-17) post-jointcal to pass the test.
        metrics : dict, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.

        Returns
        -------
        list of lsst.daf.persistence.ButlerDataRef
            The dataRefs that were processed.
        """

        result = self._runJointcalTask(nCatalogs, metrics=metrics)

        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList

        defaultBand = result.resultList[0].result.defaultBand

        def compute_statistics(refObjLoader):
            refCat = refObjLoader.loadSkyCircle(self.center, self.radius, defaultBand).refCat
            rms_result = self.jointcalStatistics.compute_rms(data_refs, refCat)
            # Make plots before testing, if requested, so we still get plots if tests fail.
            if self.do_plot:
                self._plotJointcalTask(data_refs, oldWcsList)
            return rms_result

        # we now have different astrometry/photometry refcats, so have to
        # do these calculations separately
        if self.jointcalStatistics.do_astrometry:
            refObjLoader = result.resultList[0].result.astrometryRefObjLoader
            # preserve do_photometry for the next `if`
            temp = copy.copy(self.jointcalStatistics.do_photometry)
            self.jointcalStatistics.do_photometry = False
            rms_result = compute_statistics(refObjLoader)
            self.jointcalStatistics.do_photometry = temp  # restore do_photometry

            if dist_rms_relative is not None and dist_rms_absolute is not None:
                self.assertLess(rms_result.dist_relative, dist_rms_relative)
                self.assertLess(rms_result.dist_absolute, dist_rms_absolute)

        if self.jointcalStatistics.do_photometry:
            refObjLoader = result.resultList[0].result.photometryRefObjLoader
            self.jointcalStatistics.do_astrometry = False
            rms_result = compute_statistics(refObjLoader)

            if pa1 is not None:
                self.assertLess(rms_result.pa1, pa1)

        return data_refs

    def _runJointcalTask(self, nCatalogs, metrics=None):
        """
        Run jointcalTask on nCatalogs, with the most basic tests.
        Tests for non-empty result list, and that the basic metrics are correct.

        Parameters
        ----------
        nCatalogs : int
            Number of catalogs to test on.
        metrics : dict, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.

        Returns
        -------
        pipe.base.Struct
            The structure returned by jointcalTask.run()
        """
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        if self.log_level is not None:
            self.other_args.extend(['--loglevel', 'jointcal=%s'%self.log_level])

        #  Place default configfile first so that specific subclass configfiles are applied after
        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/config.py')
        self.configfiles = [test_config] + self.configfiles

        args = [self.input_dir, '--output', self.output_dir,
                '--clobber-versions', '--clobber-config',
                '--doraise', '--configfile', *self.configfiles,
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
        data_refs : list of lsst.daf.persistence.ButlerDataRef
            The dataRefs that were processed.
        oldWcsList : list of lsst.afw.image.Wcs
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
        result : dict
            Result metric dictionary from jointcal.py
        expect : dict
            Expected metric dictionary; set a value to None to not test it.
        """
        for key in result:
            if expect[key.metric] is not None:
                value = result[key].quantity.value
                if isinstance(value, float):
                    self.assertFloatsAlmostEqual(value, expect[key.metric], msg=key.metric, rtol=1e-5)
                else:
                    self.assertEqual(value, expect[key.metric], msg=key.metric)

    def _importRepository(self, instrument, exportPath, exportFile):
        """Import a gen3 test repository into self.testDir

        Parameters
        ----------
        instrument : `str`
            Full string name for the instrument.
        exportPath : `str`
            Path to location of repository to export.
            This path must contain an `exports.yaml` file containing the
            description of the exported gen3 repo that will be imported.
        exportFile : `str`
            Filename of export data.
        """
        self.repo = os.path.join(self.output_dir, 'testrepo')

        # Make the repo and retrieve a writeable Butler
        _ = lsst.daf.butler.Butler.makeRepo(self.repo)
        butler = lsst.daf.butler.Butler(self.repo, writeable=True)
        # Register the instrument
        instrInstance = lsst.obs.base.utils.getInstrument(instrument)
        instrInstance.register(butler.registry)
        # Import the exportFile
        butler.import_(directory=exportPath, filename=exportFile,
                       transfer='symlink',
                       skip_dimensions={'instrument', 'detector', 'physical_filter'})

    def _runPipeline(self, repo, queryString=None,
                     inputCollections=None, outputCollection=None,
                     configFiles=None, configOptions=None,
                     registerDatasetTypes=False):
        """Run a pipeline via pipetask.

        Parameters
        ----------
        repo : `str`
            Gen3 Butler repository to read from/write to.
        queryString : `str`, optional
            String to use for "-d" data query. For example,
            "instrument='HSC' and tract=9697 and skymap='hsc_rings_v1'"
        inputCollections : `str`, optional
            String to use for "-i" input collections (comma delimited).
            For example, "refcats,HSC/runs/tests,HSC/calib"
        outputCollection : `str`, optional
            String to use for "-o" output collection. For example,
            "HSC/testdata/jointcal"
        configFiles : `list` [`str`], optional
            List of config files to use (with "-C").
        registerDatasetTypes : bool, optional
            Set "--register-dataset-types" when running the pipeline.

        Returns
        -------
        exit_code: `int`
            The exit code of the executed `click.testing.CliRunner`

        Raises
        ------
        RuntimeError
            Raised if the CliRunner executed pipetask failed.

        Notes
        -----
        This approach uses the Click commandline interface, which is not ideal
        for tests like this. See DM-26239 for the replacement pipetask API, and
        DM-27940 for the ticket to replace the CliRunner once that API exists.
        """
        pipelineArgs = ["run",
                        "-b", repo,
                        "-t lsst.jointcal.JointcalTask"]

        if queryString is not None:
            pipelineArgs.extend(["-d", queryString])
        if inputCollections is not None:
            pipelineArgs.extend(["-i", inputCollections])
        if outputCollection is not None:
            pipelineArgs.extend(["-o", outputCollection])
        if configFiles is not None:
            for configFile in configFiles:
                pipelineArgs.extend(["-C", configFile])
        if registerDatasetTypes:
            pipelineArgs.extend(["--register-dataset-types"])

        # TODO DM-27940: replace CliRunner with pipetask API.
        runner = click.testing.CliRunner()
        results = runner.invoke(lsst.ctrl.mpexec.cli.pipetask.cli, pipelineArgs)
        if results.exception:
            raise RuntimeError("Pipeline %s failed." % (' '.join(pipelineArgs))) from results.exception
        return results.exit_code

    def _runGen3Jointcal(self, instrumentClass, instrumentName, queryString):
        """Create a Butler repo and run jointcal on it.


        Parameters
        ----------
        instrumentClass : `str`
            The full module name of the instrument to be registered in the
            new repo. For example, "lsst.obs.subaru.HyperSuprimeCam"
        instrumentName : `str`
            The name of the instrument as it appears in the repo collections.
            For example, "HSC".
        queryString : `str`
            The query string to be run for processing. For example,
            "instrument='HSC' and tract=9697 and skymap='hsc_rings_v1'".
        """
        self._importRepository(instrumentClass,
                               self.input_dir,
                               os.path.join(self.input_dir, "exports.yaml"))
        # TODO post-RFC-741: the names of these collections will have to change
        # once testdata_jointcal is updated to reflect the collection
        # conventions in RFC-741 (no ticket for that change yet).
        inputCollections = f"refcats,{instrumentName}/testdata,{instrumentName}/calib/unbounded"
        self._runPipeline(self.repo,
                          outputCollection=f"{instrumentName}/testdata/jointcal",
                          inputCollections=inputCollections,
                          queryString=queryString,
                          registerDatasetTypes=True)
