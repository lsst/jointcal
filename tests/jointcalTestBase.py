# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import os

import lsst.afw.geom
from lsst.meas.astrom import LoadAstrometryNetObjectsTask, LoadAstrometryNetObjectsConfig

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
                   do_plot=False):
        """
        Parameters
        ----------
        center : lsst.afw.Coord
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
        """
        self._prep_reference_loader(center, radius)
        self.jointcalStatistics = utils.JointcalStatistics(match_radius)
        self.input_dir = input_dir
        self.all_visits = all_visits
        if other_args is None:
            other_args = []
        self.other_args = other_args
        self.do_plot = do_plot
        # Signal/Noise (flux/fluxSigma) for sources to be included in the RMS cross-match.
        # 100 is a balance between good centroids and enough sources.
        self.flux_limit = 100

    def tearDown(self):
        if getattr(self, 'reference', None) is not None:
            del self.reference
        if getattr(self, 'oldWcsList', None) is not None:
            del self.oldWcsList
        if getattr(self, 'jointcalTask', None) is not None:
            del self.jointcalTask

    def _prep_reference_loader(self, center, radius):
        """
        !Setup an astrometry.net reference loader.

        @param center (afw.coord) The center of the field you're testing on.
        @param radius (afw.geom.angle) The radius to load objects around center.
        """
        refLoader = LoadAstrometryNetObjectsTask(LoadAstrometryNetObjectsConfig())
        # Make a copy of the reference catalog for in-memory contiguity.
        self.reference = refLoader.loadSkyCircle(center, radius, filterName='r').refCat.copy()

    def _testJointCalTask(self, nCatalogs, relative_error, absolute_error):
        """Test parseAndRun for jointcal on nCatalogs, requiring less than some error (arcsec)."""

        visit_list = self.all_visits[:nCatalogs]
        visits = '^'.join(str(v) for v in visit_list)
        import inspect
        # the calling method is one step back on the stack: use it to specify the output repo.
        caller = inspect.stack()[1][3]  # NOTE: could be inspect.stack()[1].function in py3.5
        output_dir = os.path.join('.test', self.__class__.__name__, caller)
        args = [self.input_dir, '--output', output_dir,
                '--clobber-versions', '--clobber-config',
                '--doraise',
                '--id', 'visit=%s'%visits]
        args.extend(self.other_args)
        result = jointcal.JointcalTask.parseAndRun(args=args, doReturnResults=True)
        self.assertNotEqual(result.resultList, [], 'resultList should not be empty')
        self.dataRefs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList

        rms_rel, rms_abs = self.jointcalStatistics.compute_rms(self.dataRefs, visit_list, self.reference)
        self.assertLess(rms_rel, relative_error)
        self.assertLess(rms_abs, absolute_error)

        if self.do_plot:
            name = self.id.strip('__main__.')
            self.jointcalStatistics.make_plots(self.dataRefs, self.visitCatalogs, oldWcsList, name=name)
