# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function
from builtins import str
from builtins import object

import os
import inspect

import lsst.afw.geom
import lsst.afw.image.utils
from lsst.meas.extensions.astrometryNet import LoadAstrometryNetObjectsTask, LoadAstrometryNetObjectsConfig

from lsst.jointcal import jointcal, utils


class JointcalTestBase(object):
    """
    Base class for jointcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """

    def setUp_base(self, center, radius,
                   match_radius=0.1*lsst.afw.geom.arcseconds,
                   input_dir="",
                   all_visits=None,
                   other_args=None,
                   do_plot=False,
                   log_level=None):
        """
        Call from your child classes's setUp() to get the necessary variables built.

        Parameters
        ----------
        center : lsst.afw.IcrsCoord
            Center of the reference catalog.
        radius : lsst.afw.Angle
            Radius from center to load reference catalog objects inside.
        match_radius : lsst.afw.Angle
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
        self._prep_reference_loader(center, radius)
        self.jointcalStatistics = utils.JointcalStatistics(match_radius, verbose=True)
        self.input_dir = input_dir
        self.all_visits = all_visits
        if other_args is None:
            other_args = []
        self.other_args = other_args
        self.do_plot = do_plot
        self.log_level = log_level
        # Signal/Noise (flux/fluxSigma) for sources to be included in the RMS cross-match.
        # 100 is a balance between good centroids and enough sources.
        self.flux_limit = 100

        # Individual tests may want to tweak the config that is passed to parseAndRun().
        self.config = None

        # Append `msg` arguments to assert failures.
        self.longMessage = True

        # Ensure that the filter list is reset for each test so that we avoid
        # confusion or contamination from other instruments.
        lsst.afw.image.utils.resetFilters()

    def tearDown(self):
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

    def _prep_reference_loader(self, center, radius):
        """
        Setup an astrometry.net reference loader.

        Parameters
        ----------

        center : afw.coord
            The center of the field you're testing on.
        radius : afw.geom.angle
            The radius to load objects around center.
        """
        refLoader = LoadAstrometryNetObjectsTask(LoadAstrometryNetObjectsConfig())
        # Make a copy of the reference catalog for in-memory contiguity.
        self.reference = refLoader.loadSkyCircle(center, radius, filterName='r').refCat.copy()

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

        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5

        result = self._runJointcalTask(nCatalogs, caller, metrics=metrics)

        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)
        # Make plots before testing, if requested, so we still get plots if tests fail.
        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        if dist_rms_relative is not None and dist_rms_absolute is not None:
            self.assertLess(rms_result.dist_relative, dist_rms_relative)
            self.assertLess(rms_result.dist_absolute, dist_rms_absolute)
        if pa1 is not None:
            self.assertLess(rms_result.pa1, pa1)

        return data_refs

    def _runJointcalTask(self, nCatalogs, caller, metrics=None):
        """
        Run jointcalTask on nCatalogs, with the most basic tests.
        Tests for non-empty result list, and that the basic metrics are correct.

        Parameters
        ----------
        nCatalogs : int
            Number of catalogs to test on.
        caller : str
            Name of the calling function (to determine output directory).
        metrics : dict, optional
            Dictionary of 'metricName': value to test jointcal's result.metrics
            against.

        Returns
        -------
        pipe.base.Struct
            The structure returned by jointcalTask.run()
        """
        visits = '^'.join(str(v) for v in self.all_visits[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        if self.log_level is not None:
            self.other_args.extend(['--loglevel', 'jointcal=%s'%self.log_level])
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config',
                '--doraise',
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True, config=self.config)
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.assertEqual(result.resultList[0].exitStatus, 0)
        job = result.resultList[0].result.job
        self._test_metrics(job.measurements, metrics)

        return result

    def _plotJointcalTask(self, data_refs, oldWcsList, caller):
        """
        Plot the results of a jointcal run.

        Parameters
        ----------
        data_refs : list of lsst.daf.persistence.ButlerDataRef
            The dataRefs that were processed.
        oldWcsList : list of lsst.afw.image.Wcs
            The original WCS from each dataRef.
        caller : str
            Name of the calling function (to determine output directory).
        """
        plot_dir = os.path.join('.test', self.__class__.__name__, 'plots')
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)
        self.jointcalStatistics.make_plots(data_refs, oldWcsList, name=caller, outdir=plot_dir)
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
                if type(value) == float:
                    self.assertFloatsAlmostEqual(value, expect[key.metric], msg=key.metric, rtol=1e-5)
                else:
                    self.assertEqual(value, expect[key.metric], msg=key.metric)
