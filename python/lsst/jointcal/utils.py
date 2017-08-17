# See COPYRIGHT file at the top of the source tree.
"""
Statistics of jointcal vs. single-frame procesing and diagnostic plots.

NOTE: some of the algorithms and data structures in this code are temporary
kludges and will no longer be necessary once the following are available:
 * a composite data structure that contains all ccds from a single visit
 * an n-way matching system that preserves the separations between sources
"""
from __future__ import division, print_function, absolute_import
from builtins import zip
from builtins import object
import collections
import os

import numpy as np
from astropy import units as u

import lsst.log
import lsst.afw.table
from lsst.afw.image import fluxFromABMag, abMagFromFlux, bboxFromMetadata
from lsst.afw.geom import arcseconds

MatchDict = collections.namedtuple('MatchDict', ['relative', 'absolute'])


class JointcalStatistics(object):
    """
    Compute statistics on jointcal-processed data, and optionally generate plots.

    Notes
    -----
    Instantiate JointcalStatistics and call compute_rms() to get the relevant
    statistics for e.g. unittests, and call make_plots() to generate a suite of
    diagnostic plots.
    """

    def __init__(self, match_radius=0.1*arcseconds, flux_limit=100.0,
                 do_photometry=True, do_astrometry=True,
                 verbose=False):
        """
        Parameters
        ----------
        match_radius : lsst.afw.Angle
            match sources within this radius for RMS statistics
        flux_limit : float
            Signal/Noise (flux/fluxSigma) for sources to be included in the RMS cross-match.
            100 is a balance between good centroids and enough sources.
        do_photometry : bool, optional
            Perform calculations/make plots for photometric metrics.
        do_astrometry : bool, optional
            Perform calculations/make plots for astrometric metrics.
        verbose : bool, optional
            Print extra things
        """
        self.match_radius = match_radius
        self.flux_limit = flux_limit
        self.do_photometry = do_photometry
        self.do_astrometry = do_astrometry
        self.verbose = verbose
        self.log = lsst.log.Log.getLogger('JointcalStatistics')

    def compute_rms(self, data_refs, reference):
        """
        Match all data_refs to compute the RMS, for all detections above self.flux_limit.

        Parameters
        ----------
        data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
            A list of data refs to do the calculations between.
        reference : lsst reference catalog
            reference catalog to do absolute matching against.

        Return
        ------
        namedtuple:
            astropy.Quantity
                Post-jointcal relative RMS of the matched sources.
            astropy.Quantity
                Post-jointcal absolute RMS of matched sources.
            float
                Post-jointcal photometric repeatability (PA1 from the SRD).
        """

        # DECAM doesn't have "filter" in its registry, so we have to get filter names from VisitInfo.
        self.filters = [ref.get('calexp').getInfo().getFilter().getName() for ref in data_refs]
        self.visits_per_dataRef = [ref.dataId['visit'] for ref in data_refs]

        def compute(catalogs, calibs):
            """Compute the relative and absolute matches in distance and flux."""
            visit_catalogs = self._make_visit_catalogs(catalogs, self.visits_per_dataRef)
            catalogs = [visit_catalogs[x] for x in self.visits_per_dataRef]
            # use the first catalog as the relative reference catalog
            # NOTE: The "first" catalog depends on the original ordering of the data_refs.
            # NOTE: Thus, because I'm doing a many-1 match in _make_match_dict,
            # the number of matches (and thus the details of the match statistics)
            # will change if the data_refs are ordered differently.
            # All the more reason to use a proper n-way matcher here. See: DM-8664
            refcat = catalogs[0]
            refcalib = calibs[0]
            dist_rel, flux_rel, ref_flux_rel, source_rel = self._make_match_dict(refcat,
                                                                                 catalogs[1:],
                                                                                 calibs[1:],
                                                                                 refcalib=refcalib)
            dist_abs, flux_abs, ref_flux_abs, source_abs = self._make_match_dict(reference, catalogs, calibs)
            dist = MatchDict(dist_rel, dist_abs)
            flux = MatchDict(flux_rel, flux_abs)
            ref_flux = MatchDict(ref_flux_rel, ref_flux_abs)
            source = MatchDict(source_rel, source_abs)
            return dist, flux, ref_flux, source

        old_cats = [ref.get('src') for ref in data_refs]
        old_calibs = [ref.get('calexp').getCalib() for ref in data_refs]
        self.old_dist, self.old_flux, self.old_ref_flux, self.old_source = compute(old_cats, old_calibs)

        # Update coordinates with the new wcs, and get the new Calibs.
        new_cats = [ref.get('src') for ref in data_refs]
        new_wcss = [ref.get('wcs') for ref in data_refs]
        new_calibs = [wcs.getCalib() for wcs in new_wcss]
        if self.do_astrometry:
            for wcs, cat in zip(new_wcss, new_cats):
                # update in-place the object coordinates based on the new wcs
                lsst.afw.table.utils.updateSourceCoords(wcs.getWcs(), cat)

        self.new_dist, self.new_flux, self.new_ref_flux, self.new_source = compute(new_cats, new_calibs)

        if self.verbose:
            print('old, new relative distance matches:',
                  len(self.old_dist.relative), len(self.new_dist.relative))
            print('old, new absolute distance matches:',
                  len(self.old_dist.absolute), len(self.new_dist.absolute))
            print('old, new relative flux matches:',
                  len(self.old_flux.relative), len(self.new_flux.relative))
            print('old, new absolute flux matches:',
                  len(self.old_flux.absolute), len(self.new_flux.absolute))

        if self.do_photometry:
            self._photometric_rms()
            if self.verbose:
                print('"photometric factor" for each data ref:')
                for ref, old, new in zip(data_refs, old_calibs, new_calibs):
                    print(tuple(ref.dataId.values()), new.getFluxMag0()[0]/old.getFluxMag0()[0])
        else:
            self.new_PA1 = None

        def rms_total(data):
            """Compute the total rms across all sources."""
            total = sum(sum(dd**2) for dd in data.values())
            n = sum(len(dd) for dd in data.values())
            return np.sqrt(total/n)

        if self.do_astrometry:
            self.old_dist_total = MatchDict(*(tuple(map(rms_total, self.old_dist))*u.radian).to(u.arcsecond))
            self.new_dist_total = MatchDict(*(tuple(map(rms_total, self.new_dist))*u.radian).to(u.arcsecond))
        else:
            self.old_dist_total = MatchDict(None, None)
            self.new_dist_total = MatchDict(None, None)

        Rms_result = collections.namedtuple("Rms_result", ["dist_relative", "dist_absolute", "pa1"])
        return Rms_result(self.new_dist_total.relative, self.new_dist_total.absolute, self.new_PA1)

    def make_plots(self, data_refs, old_wcs_list,
                   name='', interactive=False, per_ccd_plot=False, outdir='.plots'):
        """
        Make plots of various quantites to help with debugging.
        Requires that `compute_rms()` was run first.

        Parameters
        ----------
        data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
            A list of data refs to do the calculations between.
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

        self.log.info("N data_refs: %d", len(data_refs))

        if self.do_photometry:
            plot_flux_distributions(plt, self.old_mag, self.new_mag,
                                    self.old_weighted_rms, self.new_weighted_rms,
                                    self.faint, self.bright, self.old_PA1, self.new_PA1,
                                    name=name, outdir=outdir)

        def rms_per_source(data):
            """Each element of data must already be the "delta" of whatever measurement."""
            return (np.sqrt([np.mean(dd**2) for dd in data.values()])*u.radian).to(u.arcsecond)

        if self.do_astrometry:
            old_dist_rms = MatchDict(*(tuple(map(rms_per_source, self.old_dist))))
            new_dist_rms = MatchDict(*(tuple(map(rms_per_source, self.new_dist))))

            self.log.info("relative RMS (old, new): {:.2e} {:.2e}".format(self.old_dist_total.relative,
                                                                          self.new_dist_total.relative))
            self.log.info("absolute RMS (old, new): {:.2e} {:.2e}".format(self.old_dist_total.absolute,
                                                                          self.new_dist_total.absolute))
            plot_rms_histogram(plt, old_dist_rms.relative, old_dist_rms.absolute,
                               new_dist_rms.relative, new_dist_rms.absolute,
                               self.old_dist_total.relative, self.old_dist_total.absolute,
                               self.new_dist_total.relative, self.new_dist_total.absolute,
                               name=name, outdir=outdir)

            plot_all_wcs_deltas(plt, data_refs, self.visits_per_dataRef, old_wcs_list,
                                per_ccd_plot=per_ccd_plot,
                                name=name, outdir=outdir)

        if interactive:
            plt.show()
            import pdb
            pdb.set_trace()

    def _photometric_rms(self, sn_cut=300, magnitude_range=3):
        """
        Compute the photometric RMS and the photometric repeatablity values (PA1).

        Parameters
        ----------
        sn_cut : float
            The minimum signal/noise for sources to be included in the PA1 calculation.
        magnitude_range : float
            The range of magnitudes above sn_cut to include in the PA1 calculation.
        """
        def rms(flux, ref_flux):
            return np.sqrt([np.mean((ref_flux[dd] - flux[dd])**2) for dd in flux])

        self.old_rms = MatchDict(*map(rms, self.old_flux, self.old_ref_flux))
        self.new_rms = MatchDict(*map(rms, self.new_flux, self.new_ref_flux))

        # we want to use the absolute fluxes for all of these calculations.
        self.old_ref = np.fromiter(self.old_ref_flux.absolute.values(), dtype=float)
        self.new_ref = np.fromiter(self.new_ref_flux.absolute.values(), dtype=float)
        self.old_mag = np.fromiter((abMagFromFlux(r) for r in self.old_ref), dtype=float)
        self.new_mag = np.fromiter((abMagFromFlux(r) for r in self.new_ref), dtype=float)

        def signal_to_noise(sources, flux_key='slot_PsfFlux_flux', sigma_key='slot_PsfFlux_fluxSigma'):
            """Compute the mean signal/noise per source from a MatchDict of SourceRecords."""
            result = np.empty(len(sources))
            for i, src in enumerate(sources.values()):
                result[i] = np.mean([x[flux_key]/x[sigma_key] for x in src])
            return result

        old_sn = signal_to_noise(self.old_source.absolute)
        # Find the faint/bright magnitude limits that are the "flat" part of the rms/magnitude relation.
        self.faint = self.old_mag[old_sn > sn_cut].max()
        self.bright = self.faint - magnitude_range
        if self.verbose:
            print("PA1 Magnitude range: {:.3f}, {:.3f}".format(self.bright, self.faint))
        old_good = (self.old_mag < self.faint) & (self.old_mag > self.bright)
        new_good = (self.new_mag < self.faint) & (self.new_mag > self.bright)
        self.old_weighted_rms = self.old_rms.absolute/self.old_ref
        self.new_weighted_rms = self.new_rms.absolute/self.new_ref
        self.old_PA1 = np.median(self.old_weighted_rms[old_good])
        self.new_PA1 = np.median(self.new_weighted_rms[new_good])

    def _make_match_dict(self, reference, visit_catalogs, calibs, refcalib=None):
        """
        Return several dicts of sourceID:[values] over the catalogs, to be used in RMS calculations.

        Parameters
        ----------
        reference : lsst.afw.table.SourceCatalog
            Catalog to do the matching against.
        visit_catalogs : list of lsst.afw.table.SourceCatalog
            Visit source catalogs (values() produced by _make_visit_catalogs)
            to cross-match against reference.
        calibs : list of lsst.afw.image.Calib
            Exposure calibs, 1-1 coorespondent with visit_catalogs.
        refcalib : lsst.afw.image.Calib or None
            Pass a Calib here to use it to compute Janskys from the reference catalog ADU slot_flux.

        Returns
        -------
        distances : dict
            dict of sourceID: array(separation distances for that source)
        fluxes : dict
            dict of sourceID: array(fluxes (Jy) for that source)
        ref_fluxes : dict
            dict of sourceID: flux (Jy) of the reference object
        sources : dict
            dict of sourceID: list(each SourceRecord that was position-matched to this sourceID)
        """

        distances = collections.defaultdict(list)
        fluxes = collections.defaultdict(list)
        ref_fluxes = {}
        sources = collections.defaultdict(list)
        if 'slot_CalibFlux_flux' in reference.schema:
            ref_flux_key = 'slot_CalibFlux_flux'
        else:
            ref_flux_key = '{}_flux'

        def get_fluxes(match):
            """Return (flux, ref_flux) or None if either is invalid."""
            # NOTE: Protect against negative fluxes: ignore this match if we find one.
            flux = match[1]['slot_CalibFlux_flux']
            if flux < 0:
                return None
            else:
                # convert to magnitudes and then Janskys, for a useable flux.
                flux = fluxFromABMag(calib.getMagnitude(flux))

            # NOTE: Have to protect against negative reference fluxes too.
            if 'slot' in ref_flux_key:
                ref_flux = match[0][ref_flux_key]
                if ref_flux < 0:
                    return None
                else:
                    ref_flux = fluxFromABMag(refcalib.getMagnitude(ref_flux))
            else:
                # a.net fluxes are already in Janskys.
                ref_flux = match[0][ref_flux_key.format(filt)]
                if ref_flux < 0:
                    return None

            Flux = collections.namedtuple('Flux', ('flux', 'ref_flux'))
            return Flux(flux, ref_flux)

        for cat, calib, filt in zip(visit_catalogs, calibs, self.filters):
            good = (cat.get('base_PsfFlux_flux')/cat.get('base_PsfFlux_fluxSigma')) > self.flux_limit
            # things the classifier called sources are not extended.
            good &= (cat.get('base_ClassificationExtendedness_value') == 0)
            matches = lsst.afw.table.matchRaDec(reference, cat[good], self.match_radius)
            for m in matches:
                if self.do_photometry:
                    flux = get_fluxes(m)
                    if flux is None:
                        continue
                    else:
                        fluxes[m[0].getId()].append(flux.flux)
                        # we can just use assignment here, since the value is always the same.
                        ref_fluxes[m[0].getId()] = flux.ref_flux

                if self.do_astrometry:
                    # Just use the computed separation distance directly.
                    distances[m[0].getId()].append(m[2])

                sources[m[0].getId()].append(m[1])
        # Convert to numpy array for easier math
        for source in distances:
            distances[source] = np.array(distances[source])
        for source in fluxes:
            fluxes[source] = np.array(fluxes[source])

        return distances, fluxes, ref_fluxes, sources

    def _make_visit_catalogs(self, catalogs, visits):
        """
        Merge all catalogs from the each visit.
        NOTE: creating this structure is somewhat slow, and will be unnecessary
        once a full-visit composite dataset is available.

        Parameters
        ----------
        catalogs : list of lsst.afw.table.SourceCatalog
            Catalogs to combine into per-visit catalogs.
        visits : list of visit id (usually int)
            list of visit identifiers, one-to-one correspondent with catalogs.

        Returns
        -------
        dict
            dict of visit: catalog of all sources from all CCDs of that visit.
        """
        visit_dict = {v: lsst.afw.table.SourceCatalog(catalogs[0].schema) for v in visits}
        for v, cat in zip(visits, catalogs):
            visit_dict[v].extend(cat)
        # We want catalog contiguity to do object selection later.
        for v in visit_dict:
            visit_dict[v] = visit_dict[v].copy(deep=True)

        return visit_dict


def plot_flux_distributions(plt, old_mag, new_mag, old_weighted_rms, new_weighted_rms,
                            faint, bright, old_PA1, new_PA1,
                            name='', outdir='.plots'):
    """Plot various distributions of fluxes and magnitudes.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with
    old_mag : np.array
        old magnitudes
    new_mag : np.array
        new magnitudes
    old_weighted_rms : np.array
        old rms weighted by the mean (rms(data)/mean(data))
    new_weighted_rms : np.array
        old rms weighted by the mean (rms(data)/mean(data))
    faint : float
        Faint end of range that PA1 was computed from.
    bright : float
        Bright end of range that PA1 was computed from.
    old_PA1 : float
        Old value of PA1, to plot as horizontal line.
    new_PA1 : float
        New value of PA1, to plot as horizontal line.
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """

    import seaborn
    seaborn.set_style('whitegrid')
    import scipy.stats

    old_color = 'blue'
    new_color = 'red'
    plt.figure()
    plt.plot(old_mag, old_weighted_rms, '.', color=old_color, label='old')
    plt.plot(new_mag, new_weighted_rms, '.', color=new_color, label='new')
    plt.axvline(faint, ls=':', color=old_color)
    plt.axvline(bright, ls=':', color=old_color)
    plt.axhline(old_PA1, ls='--', color=old_color)
    plt.axhline(new_PA1, ls='--', color=new_color)
    plt.legend(loc='upper left')
    plt.title('Where is the systematic flux rms limit?')
    plt.xlabel('magnitude')
    plt.ylabel('rms/mean per source')
    filename = os.path.join(outdir, '{}-photometry-PA1.pdf')
    plt.savefig(filename.format(name))

    plt.figure()
    seaborn.distplot(old_weighted_rms, fit=scipy.stats.lognorm, kde=False, label="old", color=old_color)
    seaborn.distplot(new_weighted_rms, fit=scipy.stats.lognorm, kde=False, label="new", color=new_color)
    plt.title('Source RMS pre/post-jointcal')
    plt.xlabel('rms(flux)/mean(flux)')
    plt.ylabel('number')
    plt.legend(loc='upper right')
    filename = os.path.join(outdir, '{}-photometry-rms.pdf')
    plt.savefig(filename.format(name))


def plot_all_wcs_deltas(plt, data_refs, visits, old_wcs_list, per_ccd_plot=False,
                        name='', outdir='.plots'):
    """
    Various plots of the difference between old and new Wcs.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with.
    data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
        A list of data refs to plot.
    visits : list of visit id (usually int)
        list of visit identifiers, one-to-one correspondent with catalogs.
    old_wcs_list : list of lsst.afw.image.wcs.Wcs
        A list of the old (pre-jointcal) WCSs, one-to-one corresponding to data_refs.
    per_ccd_plot : bool, optional
        Make per-ccd plots of the "wcs different" (warning: slow!)
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """

    plot_wcs_magnitude(plt, data_refs, visits, old_wcs_list, name, outdir=outdir)
    plot_all_wcs_quivers(plt, data_refs, visits, old_wcs_list, name, outdir=outdir)

    if per_ccd_plot:
        for i, ref in enumerate(data_refs):
            md = ref.get('calexp_md')
            dims = bboxFromMetadata(md).getDimensions()
            plot_wcs(plt, old_wcs_list[i], ref.get('wcs').getWcs(),
                     dims.getWidth(), dims.getHeight(),
                     center=(md.get('CRVAL1'), md.get('CRVAL2')), name='dataRef %d'%i,
                     outdir=outdir)


def make_xy_wcs_grid(x_dim, y_dim, wcs1, wcs2, num=50):
    """Return num x/y grid coordinates for wcs1 and wcs2."""
    x = np.linspace(0, x_dim, num)
    y = np.linspace(0, y_dim, num)
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


def plot_all_wcs_quivers(plt, data_refs, visits, old_wcs_list, name, outdir='.plots'):
    """
    Make quiver plots of the WCS deltas for each CCD in each visit.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with.
    data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
        A list of data refs to plot.
    visits : list of visit id (usually int)
        list of visit identifiers, one-to-one correspondent with catalogs.
    old_wcs_list : list of lsst.afw.image.wcs.Wcs
        A list of the old (pre-jointcal) WCSs, one-to-one corresponding to data_refs.
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """

    for visit in visits:
        fig = plt.figure()
        # fig.set_tight_layout(True)
        ax = fig.add_subplot(111)
        for old_wcs, ref in zip(old_wcs_list, data_refs):
            if ref.dataId['visit'] != visit:
                continue
            md = ref.get('calexp_md')
            dims = bboxFromMetadata(md).getDimensions()
            Q = plot_wcs_quivers(ax, old_wcs, ref.get('wcs').getWcs(),
                                 dims.getWidth(), dims.getHeight())
            # TODO: add CCD bounding boxes to plot once DM-5503 is finished.
            # TODO: add a circle for the full focal plane.
        length = (0.1*u.arcsecond).to(u.radian).value
        ax.quiverkey(Q, 0.9, 0.95, length, '0.1 arcsec', coordinates='figure', labelpos='W')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.title('visit: {}'.format(visit))
        filename = os.path.join(outdir, '{}-{}-quivers.pdf')
        plt.savefig(filename.format(name, visit))


def plot_wcs_quivers(ax, wcs1, wcs2, x_dim, y_dim):
    """
    Plot the delta between wcs1 and wcs2 as vector arrows.

    Parameters
    ----------
    ax : matplotlib.axis
        Matplotlib axis instance to plot to.
    wcs1 : lsst.afw.image.wcs.Wcs
        First WCS to compare.
    wcs2 : lsst.afw.image.wcs.Wcs
        Second WCS to compare.
    x_dim : int
        Size of array in X-coordinate to make the grid over.
    y_dim : int
        Size of array in Y-coordinate to make the grid over.
    """

    x1, y1, x2, y2 = make_xy_wcs_grid(x_dim, y_dim, wcs1, wcs2)
    uu = x2 - x1
    vv = y2 - y1
    return ax.quiver(x1, y1, uu, vv, units='x', pivot='tail', scale=1e-3, width=1e-5)


def plot_wcs_magnitude(plt, data_refs, visits, old_wcs_list, name, outdir='.plots'):
    """Plot the magnitude of the WCS change between old and new visits as a heat map.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with.
    data_refs : list of lsst.daf.persistence.butlerSubset.ButlerDataRef
        A list of data refs to plot.
    visits : list of visit id (usually int)
        list of visit identifiers, one-to-one correspondent with catalogs.
    old_wcs_list : list of lsst.afw.image.wcs.Wcs
        A list of the old (pre-jointcal) WCSs, one-to-one corresponding to data_refs.
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """
    for visit in visits:
        fig = plt.figure()
        fig.set_tight_layout(True)
        ax = fig.add_subplot(111)
        # Start min/max at the "opposite" ends so they always get the first valid value.
        xmin = np.inf
        ymin = np.inf
        xmax = -np.inf
        ymax = -np.inf
        for old_wcs, ref in zip(old_wcs_list, data_refs):
            if ref.dataId['visit'] != visit:
                continue
            md = ref.get('calexp_md')
            dims = bboxFromMetadata(md).getDimensions()
            x1, y1, x2, y2 = make_xy_wcs_grid(dims.getWidth(), dims.getHeight(),
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
        plt.title('visit: {}'.format(visit))
        filename = os.path.join(outdir, '{}-{}-heatmap.pdf')
        plt.savefig(filename.format(name, visit))


def plot_wcs(plt, wcs1, wcs2, x_dim, y_dim, center=(0, 0), name="", outdir='.plots'):
    """Plot the "distortion map": wcs1-wcs2 delta of points in the CCD grid.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with.
    wcs1 : lsst.afw.image.wcs.Wcs
        First WCS to compare.
    wcs2 : lsst.afw.image.wcs.Wcs
        Second WCS to compare.
    x_dim : int
        Size of array in X-coordinate to make the grid over.
    y_dim : int
        Size of array in Y-coordinate to make the grid over.
    center : tuple, optional
        Center of the data, in on-chip coordinates.
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """

    plt.figure()

    x1, y1, x2, y2 = make_xy_wcs_grid(x_dim, y_dim, wcs1, wcs2, num=50)
    plt.plot((x1 - x2) + center[0], (y1 - y2) + center[1], '-')
    plt.xlabel('delta RA (arcsec)')
    plt.ylabel('delta Dec (arcsec)')
    plt.title(name)
    filename = os.path.join(outdir, '{}-wcs.pdf')
    plt.savefig(filename.format(name))


def plot_rms_histogram(plt, old_rms_relative, old_rms_absolute,
                       new_rms_relative, new_rms_absolute,
                       old_rel_total, old_abs_total, new_rel_total, new_abs_total,
                       name="", outdir='.plots'):
    """Plot histograms of the source separations and their RMS values.

    Parameters
    ----------
    plt : matplotlib.pyplot instance
        pyplot instance to plot with.
    old_rms_relative : np.array
        old relative rms/star
    old_rms_absolute : np.array
        old absolute rms/star
    new_rms_relative : np.array
        new relative rms/star
    new_rms_absolute : np.array
        new absolute rms/star
    old_rel_total : float
        old relative rms over all stars
    old_abs_total : float
        old absolute rms over all stars
    new_rel_total : float
        new relative rms over all stars
    new_abs_total : float
        new absolute rms over all stars
    name : str
        Name to include in plot titles and save files.
    outdir : str, optional
        Directory to write the saved plots to.
    """
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
