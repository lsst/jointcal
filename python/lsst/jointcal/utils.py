# See COPYRIGHT file at the top of the source tree.
"""
Statistics of jointcal vs. single-frame procesing and diagnostic plots.

NOTE: some of the algorithms and data structures in this code are temporary
kludges and will no longer be necessary once the following are available:
 * a composite data structure that contains all ccds from a single visit
 * an n-way matching system that preserves the separations between sources
"""
import collections
import os

import numpy as np
from astropy import units as u

import lsst.log
import lsst.afw.table
from lsst.afw.geom import arcseconds


class JointcalStatistics(object):
    """Compute statistics on jointcal-processed data, and optionally generate plots."""

    def __init__(self, match_radius=0.1*arcseconds, flux_limit=100.):
        """
        Parameters
        ----------
        match_radius : lsst.afw.Angle
            match sources within this radius for RMS statistics
        flux_limit : float
            Signal/Noise (flux/fluxSigma) for sources to be included in the RMS cross-match.
            100 is a balance between good centroids and enough sources.
        """
        self.match_radius = match_radius
        self.flux_limit = flux_limit
        self.log = lsst.log.Log.getLogger('JointcalStatistics')

    def compute_rms(self, data_refs, visit_list, reference):
        """
        Match all data_refs to compute the RMS, for all detections above self.flux_limit.

        Parameters
        ----------
        data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
            A list of data refs to do the calculations between.
        visit_list : list of visit id (usually int)
            list of visit identifiers to do the catalog merge on.
        reference : lsst reference catalog
            reference catalog to do absolute matching against.

        Return
        ------
        geom.Angle
            New relative RMS of the matched sources.
        geom.Angle
            New absolute RMS of matched sources.
        """
        visits_per_dataRef = [dataRef.dataId['visit'] for dataRef in data_refs]

        def compute(catalogs):
            """Compute the relative and absolute RMS."""
            visit_catalogs = self._make_visit_catalogs(catalogs, visits_per_dataRef, visit_list)
            refcat = visit_catalogs.values()[0]  # use the first catalog as the relative reference catalog
            relative = self._make_match_dict(refcat, visit_catalogs.values()[1:])
            absolute = self._make_match_dict(reference, visit_catalogs.values())
            return visit_catalogs, relative, absolute

        old_cats = [dataRef.get('src') for dataRef in data_refs]
        visit_catalogs, self.old_relative, self.old_absolute = compute(old_cats)

        # Update coordinates with the new wcs including distortion corrections.
        new_cats = []
        for dataRef in data_refs:
            new_cats.append(dataRef.get('src'))
            lsst.afw.table.utils.updateSourceCoords(dataRef.get('wcs').getWcs(), new_cats[-1])
        visistCatalogs, self.new_relative, self.new_absolute = compute(new_cats)

        self.old_rel_total = rms_total(self.old_relative)
        self.new_rel_total = rms_total(self.new_relative)
        self.old_abs_total = rms_total(self.old_absolute)
        self.new_abs_total = rms_total(self.new_absolute)

        return self.new_rel_total, self.new_abs_total

    def make_plots(self, data_refs, visit_catalogs, old_wcs_list, name='',
                   interactive=False, per_ccd_plot=False, outdir='.plots'):
        """
        Make plots of various quantites to help with debugging.
        Requires that `compute_rms()` was run first.

        Parameters
        ----------
        data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
            A list of data refs to do the calculations between.
        visit_catalogs : list of lsst.afw.table.SourceCatalog
            visit source catalogs (values() produced by _make_visit_catalogs) to cross-match.
        old_wcs_list : list of lsst.afw.image.wcs.Wcs
            A list of the old (pre-jointcal) WCSs, one-to-one corresponding to data_refs.
        name : str
            Name to include in plot titles and save files.
        interactive : bool
            Turn on matplotlib interactive mode and drop into a debugger when
            plotting is finished. Otherwise, use a non-interactive backend.
        per_ccd_plot : bool
            Plot the WCS per CCD (takes longer and generates many plots for a large camera)
        outdir : str
            directory to save plots to.
        """
        import matplotlib
        if not interactive:
            # Use a non-interactive backend for faster plotting.
            matplotlib.use('pdf')

        import matplotlib.pyplot as plt
        import astropy.visualization
        # make quantities behave nicely when plotted.
        astropy.visualization.quantity_support()
        if interactive:
            plt.ion()

        old_rms_relative = rms_per_source(self.old_relative)
        old_rms_absolute = rms_per_source(self.old_absolute)
        new_rms_relative = rms_per_source(self.new_relative)
        new_rms_absolute = rms_per_source(self.new_absolute)
        self.log.info("N data_refs: %d", len(data_refs))
        self.log.info("relative RMS (old, new): {:.2e} {:.2e}".format(self.old_rel_total, self.new_rel_total))
        self.log.info("absolute RMS (old, new): {:.2e} {:.2e}".format(self.old_abs_total, self.new_abs_total))
        plot_rms_histogram(plt, old_rms_relative, old_rms_absolute,
                           new_rms_relative, new_rms_absolute,
                           self.old_rel_total, self.old_abs_total,
                           self.new_rel_total, self.new_abs_total, name, outdir=outdir)

        plot_all_wcs_deltas(plt, data_refs, visit_catalogs, old_wcs_list, name,
                            outdir=outdir, per_ccd_plot=per_ccd_plot)

        if interactive:
            plt.show()
            import pdb
            pdb.set_trace()

    def _make_match_dict(self, reference, visit_catalogs):
        """
        Return a dict of sourceID:[distances] over the catalogs, for RMS calculations.

        Parameters
        ----------
        reference : lsst.afw.table.SourceCatalog
            Catalog to do the matching against.
        visit_catalogs : list of lsst.afw.table.SourceCatalog
            Visit source catalogs (values() produced by _make_visit_catalogs)
            to cross-match against reference.

        Returns
        -------
        dict
            dict of sourceID:list(separation distances for that source)
        """

        deltas = collections.defaultdict(list)
        for cat in visit_catalogs:
            good = (cat.get('base_PsfFlux_flux')/cat.get('base_PsfFlux_fluxSigma')) > self.flux_limit
            # things the classifier called sources are not extended.
            good &= (cat.get('base_ClassificationExtendedness_value') == 0)
            matches = lsst.afw.table.matchRaDec(reference, cat[good], self.match_radius)
            for m in matches:
                deltas[m[0].getId()].append(m[2])
        return deltas

    def _make_visit_catalogs(self, catalogs, visits, visit_list):
        """
        Merge all catalogs from the each visit.
        NOTE: creating this structure is somewhat slow, and will be unnecessary
        once a full-visit composite dataset is available.

        Parameters
        ----------
        catalogs : list of lsst.afw.table.SourceCatalog
            Catalogs to combine into per-visit catalogs.
        visits : list of visit id (usually int)
            list of visit identifiers, one-to-one correspondant with catalogs.
        visit_list : list of visit id (usually int)
            list of visit identifiers to do the catalog merge on (a proper subset of visits).

        Returns
        -------
        dict
            dict of visit: catalog of all sources from all CCDs of that visit.
        """
        visit_dict = {v: lsst.afw.table.SourceCatalog(catalogs[0].schema) for v in visit_list}
        for v, cat in zip(visits, catalogs):
            visit_dict[v].extend(cat)
        # We want catalog contiguity to do object selection later.
        for v in visit_dict:
            visit_dict[v] = visit_dict[v].copy(deep=True)

        return visit_dict


def rms_per_source(data):
    """
    Compute the RMS for each catalog source.

    Parameters
    ----------
    data : dict
        dict of sourceID:list(measurement deltas for that source)

    Returns
    -------
    np.array of astropy.units.arcsecond
        the RMS for each sourceID (no guaranteed order)
    """
    return (np.sqrt(np.array([np.mean(np.array(dd)**2) for dd in data.values()]))*u.radian).to(u.arcsecond)


def rms_total(data):
    """
    Compute the total rms across all sources.
    Parameters
    ----------
    data : dict
        dict of sourceID:list(measurement deltas for that source)

    Returns
    -------
    astropy.units.arcsecond
        the total rms across all sources
    """
    total = sum(sum(np.array(dd)**2) for dd in data.values())
    n = sum(len(dd) for dd in data.values())
    return (np.sqrt(total/n) * u.radian).to(u.arcsecond)


def plot_all_wcs_deltas(plt, data_refs, visit_catalogs, old_wcs_list, name,
                        per_ccd_plot=False, outdir='.plots'):
    """Various plots of the difference between old and new Wcs."""

    plot_wcs_magnitude(plt, data_refs, visit_catalogs, old_wcs_list, name, outdir=outdir)
    plot_all_wcs_quivers(plt, data_refs, visit_catalogs, old_wcs_list, name, outdir=outdir)

    if per_ccd_plot:
        for i, ref in enumerate(data_refs):
            md = ref.get('calexp_md')
            plot_wcs(plt, old_wcs_list[i], ref.get('wcs').getWcs(),
                     dim=(md.get('NAXIS1'), md.get('NAXIS1')),
                     center=(md.get('CRVAL1'), md.get('CRVAL2')), name='dataRef %d'%i,
                     outdir=outdir)


def make_xy_wcs_grid(dim, wcs1, wcs2, num=50):
    """Return num x/y grid coordinates for wcs1 and wcs2."""
    x = np.linspace(0, dim[0], num)
    y = np.linspace(0, dim[1], num)
    x1, y1 = wcs_convert(x, y, wcs1)
    x2, y2 = wcs_convert(x, y, wcs2)
    return x1, y1, x2, y2


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


def plot_all_wcs_quivers(plt, data_refs, visit_catalogs, old_wcs_list, name, outdir='.plots'):
    """Make quiver plots of the WCS deltas for each CCD in each visit."""

    for cat in visit_catalogs:
        fig = plt.figure()
        # fig.set_tight_layout(True)
        ax = fig.add_subplot(111)
        for old_wcs, ref in zip(old_wcs_list, data_refs):
            if ref.dataId['visit'] != cat:
                continue
            md = ref.get('calexp_md')
            Q = plot_wcs_quivers(ax, old_wcs, ref.get('wcs').getWcs(),
                                 dim=(md.get('NAXIS1'), md.get('NAXIS2')))
            # TODO: add CCD bounding boxes to plot once DM-5503 is finished.
            # TODO: add a circle for the full focal plane.
        length = (0.1*u.arcsecond).to(u.radian).value
        ax.quiverkey(Q, 0.9, 0.95, length, '0.1 arcsec', coordinates='figure', labelpos='W')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.title('visit: {}'.format(cat))
        filename = os.path.join(outdir, '{}-{}-quivers.pdf')
        plt.savefig(filename.format(name, cat))


def plot_wcs_quivers(ax, wcs1, wcs2, dim):
    """Plot the delta between wcs1 and wcs2 as vector arrows."""

    x1, y1, x2, y2 = make_xy_wcs_grid(dim, wcs1, wcs2)
    uu = x2 - x1
    vv = y2 - y1
    return ax.quiver(x1, y1, uu, vv, units='x', pivot='tail', scale=1e-3, width=1e-5)


def plot_wcs_magnitude(plt, data_refs, visit_catalogs, old_wcs_list, name, outdir='.plots'):
    """Plot the magnitude of the WCS change between old and new visits as a heat map."""
    for cat in visit_catalogs:
        fig = plt.figure()
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)
        # Start min/max at the "opposite" ends so they always get the first valid value.
        xmin = np.inf
        ymin = np.inf
        xmax = -np.inf
        ymax = -np.inf
        for old_wcs, ref in zip(old_wcs_list, data_refs):
            if ref.dataId['visit'] != cat:
                continue
            md = ref.get('calexp_md')
            x1, y1, x2, y2 = make_xy_wcs_grid((md.get('NAXIS1'), md.get('NAXIS2')),
                                              old_wcs, ref.get('wcs').getWcs())
            uu = x2 - x1
            vv = y2 - y1
            extent = (x1[0, 0], x1[-1, -1], y1[0, 0], y1[-1, -1])
            xmin = min(x1.min(), xmin)
            ymin = min(y1.min(), ymin)
            xmax = max(x1.max(), xmax)
            ymax = max(y1.max(), ymax)
            magnitude = (np.linalg.norm((uu, vv), axis=0)*u.radian).to(u.arcsecond).value
            img = ax.imshow(magnitude, vmin=0, vmax=0.3,
                            aspect='auto', extent=extent, cmap=plt.get_cmap('magma'))
            # TODO: add CCD bounding boxes to the plot once DM-5503 is finished.
            # TODO: add a circle for the full focal plane.

        # We're reusing only one of the returned images here for colorbar scaling,
        # but it doesn't matter because we set vmin/vmax so they are all scaled the same.
        cbar = plt.colorbar(img)
        cbar.ax.set_ylabel('distortion (arcseconds)')
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.title('visit: {}'.format(cat))
        filename = os.path.join(outdir, '{}-{}-heatmap.pdf')
        plt.savefig(filename.format(name, cat))


def plot_wcs(plt, wcs1, wcs2, dim, center=(0, 0), name="", outdir='.plots'):
    """Plot the "distortion map": wcs1-wcs2 delta of points in the CCD grid."""

    plt.figure()

    x1, y1, x2, y2 = make_xy_wcs_grid(dim, wcs1, wcs2, num=50)
    plt.plot((x1 - x2) + center[0], (y1 - y2) + center[1], '-')
    plt.xlabel('delta RA (arcsec)')
    plt.ylabel('delta Dec (arcsec)')
    plt.title(name)
    filename = os.path.join(outdir, '{}-wcs.pdf')
    plt.savefig(filename.format(name))


def plot_rms_histogram(plt, old_rms_relative, old_rms_absolute,
                       new_rms_relative, new_rms_absolute,
                       old_rel_total, old_abs_total, new_rel_total, new_abs_total,
                       name, outdir='.plots'):
    """Plot histograms of the source separations and their RMS values."""
    plt.figure()

    color_rel = 'black'
    ls_old = 'dotted'
    color_abs = 'green'
    ls_new = 'dashed'
    plotOptions = {'lw': 2, 'range': (0, 0.1)*u.arcsecond, 'normed': True,
                   'bins': 30, 'histtype': 'step'}

    plt.title('relative vs. absolute: %d vs. %d'%(len(old_rms_relative), len(old_rms_absolute)))

    plt.hist(old_rms_absolute, color=color_abs, ls=ls_old, label='old abs', **plotOptions)
    plt.hist(new_rms_absolute, color=color_abs, ls=ls_new, label='new abs', **plotOptions)

    plt.hist(old_rms_relative, color=color_rel, ls=ls_old, label='old rel', **plotOptions)
    plt.hist(new_rms_relative, color=color_rel, ls=ls_new, label='new rel', **plotOptions)

    plt.axvline(x=old_abs_total.value, linewidth=1.5, color=color_abs, ls=ls_old)
    plt.axvline(x=new_abs_total.value, linewidth=1.5, color=color_abs, ls=ls_new)
    plt.axvline(x=old_rel_total.value, linewidth=1.5, color=color_rel, ls=ls_old)
    plt.axvline(x=new_rel_total.value, linewidth=1.5, color=color_rel, ls=ls_new)

    plt.xlim(plotOptions['range'])
    plt.xlabel('arcseconds')
    plt.legend(loc='best')
    filename = os.path.join(outdir, '{}-histogram.pdf')
    plt.savefig(filename.format(name))
