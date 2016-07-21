# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import numpy as np
import collections
import os

from lsst.meas.astrom import LoadAstrometryNetObjectsTask, LoadAstrometryNetObjectsConfig
from lsst.afw import geom, table

from lsst.jointcal import jointcal

arcsec_per_radian = 206265.


class JointcalTestBase(object):
    """
    Base class for jointcal tests, to genericize some test running and setup.

    Derive from this first, then from TestCase.
    """

    def setUp(self):
        self.doPlot = False  # don't make plots unless specifically requested
        self.matchRadius = 0.1*geom.arcseconds  # match sources within 0.1" for RMS statistics

    def tearDown(self):
        if getattr(self, 'jointcalTask', None) is not None:
            del self.jointcalTask
        if getattr(self, 'jointcalTask', None) is not None:
            del self.catalogs

    def _prep_reference_loader(self, center, radius):
        """
        Setup an astrometry.net reference loader.

        @param center (afw.coord) The center of the field you're testing on.
        @param radius (afw.geom.angle) The radius to load objects around center.
        """
        refLoader = LoadAstrometryNetObjectsTask(LoadAstrometryNetObjectsConfig())
        # Make a copy of the reference catalog for in-memory contiguity.
        self.reference = refLoader.loadSkyCircle(center, radius, filterName='r').refCat.copy()

    def _testJointCalTask_run(self, nCatalogs, relative_error, absolute_error):
        """Test jointcal.run() on nCatalogs, requiring less than some error (arcsec)."""

        result = self.jointcalTask.run(self.catalogs[:nCatalogs])
        self.dataRefs = result.dataRefs
        self.old_wcss = result.old_wcss

        rms_rel, rms_abs = self.compute_rms(self.catalogs[:nCatalogs], fluxlimit=100)
        self.assertLess(rms_rel, relative_error/arcsec_per_radian)
        self.assertLess(rms_abs, absolute_error)

    def _testJointCalTask(self, nCatalogs, relative_error, absolute_error):
        """Test parseAndRun for jointcal on nCatalogs, requiring less than some error (arcsec)."""

        visits = '^'.join(str(v) for v in self.visitList[:nCatalogs])
        output_dir = os.path.join('.test', self.__class__.__name__)
        result = jointcal.JointcalTask.parseAndRun(args=[self.input_dir, '--output', output_dir,
                                                         '--clobber-versions', '--clobber-config',
                                                         '--doraise',
                                                         '--id', 'visit=%s'%visits],
                                                   doReturnResults=True)
        self.dataRefs = result.resultList[0].result.dataRefs
        self.old_wcss = result.resultList[0].result.old_wcss

        rms_rel, rms_abs = self.compute_rms(self.dataRefs, fluxlimit=100)
        self.assertLess(rms_rel, relative_error/arcsec_per_radian)
        self.assertLess(rms_abs, absolute_error)

    def compute_rms(self, dataRefs, fluxlimit=100):
        """
        Match all dataRefs to compute the RMS, for all detections above some flux limit.

        @param dataRefs (list) the dataRefs to do the relative matching between
        @param fluxlimit (float) SN (flux/fluxSigma) for sources to be included in match.

        @return (float, float) relative, absolute RMS of stars, in radians
        """

        visit_list = [dataRef.dataId['visit'] for dataRef in dataRefs]

        old_cats = [dataRef.get('src') for dataRef in dataRefs]
        visitCatalogs = self._make_visit_catalogs(old_cats, visit_list)
        refSrc = visitCatalogs.values()[0]  # use the first catalog as the relative reference catalog
        old_relative = self._make_match_dict(refSrc, visitCatalogs.values()[1:], fluxlimit=fluxlimit)
        old_absolute = self._make_match_dict(self.reference, visitCatalogs.values(), fluxlimit=fluxlimit)

        # Update coordinates with the new wcs including distortion corrections.
        new_cats = []
        for dataRef in dataRefs:
            new_cats.append(dataRef.get('src'))
            table.utils.updateSourceCoords(dataRef.get('wcs').getWcs(), new_cats[-1])
        visitCatalogs = self._make_visit_catalogs(new_cats, visit_list)
        refSrc = visitCatalogs.values()[0]  # use the first catalog as the relative reference catalog
        new_relative = self._make_match_dict(refSrc, visitCatalogs.values()[1:], fluxlimit=fluxlimit)
        new_absolute = self._make_match_dict(self.reference, visitCatalogs.values(), fluxlimit=fluxlimit)

        old_rel_total = rms_total(old_relative)
        new_rel_total = rms_total(new_relative)
        old_abs_total = rms_total(old_absolute)
        new_abs_total = rms_total(new_absolute)

        if self.doPlot:
            self._do_plots(dataRefs, visitCatalogs,
                           old_relative, old_absolute, new_relative, new_absolute,
                           old_rel_total, old_abs_total, new_rel_total, new_abs_total)

        return new_rel_total, new_abs_total

    def _make_match_dict(self, reference, visitCatalogs, fluxlimit=100):
        """
        Return a dict of starID:[distances] over the catalogs, for RMS calculations.

        @param reference (SourceCatalog) catalog to do the matching against
        @param visits (list) visit source catalogs (from _make_visit_catalogs) to cross-match
        @param fluxlimit (float) SN (flux/fluxSigma) for sources to be included in match

        @return (dict) dict of starID:list(measurement deltas for that star)
        """

        deltas = collections.defaultdict(list)
        for cat in visitCatalogs:
            good = (cat.get('base_PsfFlux_flux')/cat.get('base_PsfFlux_fluxSigma')) > fluxlimit
            # things the classifier called stars are not extended.
            good &= (cat.get('base_ClassificationExtendedness_value') == 0)
            matches = table.matchRaDec(reference, cat[good], self.matchRadius)
            for m in matches:
                deltas[m[0].getId()].append(m[2])
        return deltas

    def _make_visit_catalogs(self, catalogs, visits):
        """
        Return a dict of visit: catalog of all sources from all CCDs of that visit.

        @param catalogs (list of SourceCatalogs) catalogs to combine into per-visit catalogs.
        @param visits (list) list of visit identifiers, one-to-one correspondant with catalogs.
        """
        visit_dict = {v: table.SourceCatalog(catalogs[0].schema) for v in self.visitList}
        for v, cat in zip(visits, catalogs):
            visit_dict[v].extend(cat)
        # We want catalog contiguity to do object selection later.
        for v in visit_dict:
            visit_dict[v] = visit_dict[v].copy(deep=True)

        return visit_dict

    def _do_plots(self, dataRefs, visitCatalogs,
                  old_relative, old_absolute, new_relative, new_absolute,
                  old_rel_total, old_abs_total, new_rel_total, new_abs_total):
        """Make plots of various quantites to help with debugging."""
        import matplotlib.pyplot as plt
        plt.ion()

        plot_all_wcs_deltas(plt, dataRefs, visitCatalogs, self.old_wcss)

        old_rms_relative = rms_per_star(old_relative)
        old_rms_absolute = rms_per_star(old_absolute)
        new_rms_relative = rms_per_star(new_relative)
        new_rms_absolute = rms_per_star(new_absolute)
        print(len(dataRefs))
        print("relative RMS (old, new):", old_rel_total*arcsec_per_radian,
              new_rel_total*arcsec_per_radian)
        print("absolute RMS (old, new):", old_abs_total*arcsec_per_radian,
              new_abs_total*arcsec_per_radian)
        plot_rms_histogram(plt, old_rms_relative, old_rms_absolute, new_rms_relative, new_rms_absolute,
                           old_rel_total, old_abs_total, new_rel_total, new_abs_total)

        # So one can muck-about with things after plotting...
        import ipdb
        ipdb.set_trace()


def rms_per_star(data):
    """
    Compute the RMS/catalog star

    @param data (dict) dict of starID:list(measurement deltas for that star)

    @return (np.array) the rms for each starID (no guaranteed order)
    """
    result = np.empty(len(data))
    for i, k in enumerate(data):
        result[i] = np.sqrt(np.mean(np.array(data[k])**2))
    return result


def rms_total(data):
    """
    Compute the total rms across all stars.

    @param data (dict) dict of starID:list(measurement deltas for that star)

    @return (int) the total rms across all stars
    """
    n = 0
    total = 0
    for i, k in enumerate(data):
        total += sum(np.array(data[k])**2)
        n += len(data[k])
    return np.sqrt(total/n)


def plot_all_wcs_deltas(plt, dataRefs, visitCatalogs, old_wcss):
    """Various plots of the difference between old and new Wcs."""

    plot_all_wcs_quivers(plt, dataRefs, visitCatalogs, old_wcss)

    # NOTE: uncomment this block to get per-ccd plots of the "distortion map".
    # for i, ref in enumerate(dataRefs):
    #     md = ref.get('calexp_md')
    #     plot_wcs(plt, old_wcss[i], ref.get('wcs').getWcs(),
    #              dim=(md.get('NAXIS1'), md.get('NAXIS1')),
    #              center=(md.get('CRVAL1'), md.get('CRVAL2')), name='dataRef %d'%i)


def wcs_convert(xv, yv, wcs):
    """Convert two arrays of x/y points into an on-sky grid."""
    xout = np.zeros((xv.shape[0], yv.shape[0]))
    yout = np.zeros((xv.shape[0], yv.shape[0]))
    for i, x in enumerate(xv):
        for j, y in enumerate(yv):
            sky = wcs.pixelToSky(x, y).toFk5()
            xout[i, j] = sky.getRa()
            yout[i, j] = sky.getDec()
    return xout, yout


def plot_all_wcs_quivers(plt, dataRefs, visitCatalogs, old_wcss):
    """Make quiver plots of the WCS deltas for each CCD in each visit."""

    for cat in visitCatalogs:
        fig = plt.figure()
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)
        for old_wcs, ref in zip(old_wcss, dataRefs):
            if ref.dataId['visit'] == cat:
                md = ref.get('calexp_md')
                Q = plot_wcs_quivers(ax, old_wcs, ref.get('wcs').getWcs(),
                                     dim=(md.get('NAXIS1'), md.get('NAXIS2')))
                # TODO: add CCD bounding boxes to plot once DM-5503 is finished.
                # TODO: add a circle for the full focal plane.
        ax.quiverkey(Q, 0.9, 0.95, 0.1/arcsec_per_radian, '0.1 arcsec', coordinates='figure')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.title('visit: {}'.format(cat))


def plot_wcs_quivers(ax, wcs1, wcs2, dim=(1000, 1000)):
    """Plot the delta between wcs1 and wcs2 as vector arrows."""

    x = np.linspace(0, dim[0], 50)
    y = np.linspace(0, dim[0], 50)
    x1, y1 = wcs_convert(x, y, wcs1)
    x2, y2 = wcs_convert(x, y, wcs2)
    u = x2-x1
    v = y2-y1
    Q = ax.quiver(x1, y1, u, v, units='x', pivot='tail', scale=1e-3, width=5e-6)
    return Q


def plot_wcs(plt, wcs1, wcs2, dim=(1000, 1000), center=(0, 0), name=""):
    """Plot the "distortion map": wcs1-wcs2 delta of points in the CCD grid."""

    plt.figure()

    x = np.linspace(0, dim[0], 50)
    y = np.linspace(0, dim[0], 50)
    x1, y1 = wcs_convert(x, y, wcs1)
    x2, y2 = wcs_convert(x, y, wcs2)

    plt.plot((x1-x2)*arcsec_per_radian + center[0], (y1-y2)*arcsec_per_radian + center[1], '-')
    plt.xlabel('delta RA (arcsec)')
    plt.ylabel('delta Dec (arcsec)')
    plt.title(name)


def plot_rms_histogram(plt, old_rms_relative, old_rms_absolute,
                       new_rms_relative, new_rms_absolute,
                       old_rel_total, old_abs_total, new_rel_total, new_abs_total):
    """Plot histograms of the star separations and their RMS values."""
    plt.ion()
    plt.figure()

    color_rel = 'black'
    ls_old = 'dotted'
    color_abs = 'blue'
    ls_new = 'dashed'
    lw = 2

    plt.title('relative vs. absolute: %d vs. %d'%(len(old_rms_relative), len(old_rms_absolute)))

    plt.hist(old_rms_absolute*arcsec_per_radian, color=color_abs, ls=ls_old, label='old abs',
             histtype='step', bins=30, lw=lw)
    plt.hist(new_rms_absolute*arcsec_per_radian, color=color_abs, ls=ls_new, label='new abs',
             histtype='step', bins=30, lw=lw)

    plt.hist(old_rms_relative*arcsec_per_radian, color=color_rel, ls=ls_old, label='old rel',
             histtype='step', bins=30, lw=lw)
    plt.hist(new_rms_relative*arcsec_per_radian, color=color_rel, ls=ls_new, label='new rel',
             histtype='step', bins=30, lw=lw)

    plt.axvline(x=old_abs_total*arcsec_per_radian, linewidth=1.5, color=color_abs, ls=ls_old)
    plt.axvline(x=new_abs_total*arcsec_per_radian, linewidth=1.5, color=color_abs, ls=ls_new)
    plt.axvline(x=old_rel_total*arcsec_per_radian, linewidth=1.5, color=color_rel, ls=ls_old)
    plt.axvline(x=new_rel_total*arcsec_per_radian, linewidth=1.5, color=color_rel, ls=ls_new)

    plt.xlabel('arcseconds')
    plt.legend(loc='best')
