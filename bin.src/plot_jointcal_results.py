#!/usr/bin/env python
# See COPYRIGHT file at the top of the source tree.
"""
Generate plots from the results of a jointcal run.

Requires a butler repository that has had jointcal successfully run on it, which
produces a new "wcs" component of the repository.
"""

import os

import lsst.daf.persistence
from lsst.meas.astrom import LoadAstrometryNetObjectsTask, LoadAstrometryNetObjectsConfig
import lsst.afw.image
from lsst.afw.geom import degrees

from lsst.jointcal import utils


def get_valid_dataIds(butler, dataset_type='wcs'):
    """Return a list of all dataIds that exist in this butler."""
    data_ids = []
    for data_ref in butler.subset(dataset_type):
        if data_ref.datasetExists():
            data_ids.append(data_ref.dataId)
    return data_ids


def prep_reference_loader(center, radius):
    """
    Return an astrometry.net reference loader.

    Parameters
    ----------
    center: afw.coord
        The center of the field you're testing on.
    radius: afw.geom.angle
        The radius to load objects around center.
    """
    refLoader = LoadAstrometryNetObjectsTask(LoadAstrometryNetObjectsConfig())
    # Make a copy of the reference catalog for in-memory contiguity.
    return refLoader.loadSkyCircle(center, radius, filterName='r').refCat.copy()


def get_old_wcs_list(data_refs):
    """Return a list of the original WCSs of the data_refs."""
    result = []
    for data_ref in data_refs:
        md = data_ref.get("calexp_md", immediate=True)
        tanwcs = lsst.afw.image.TanWcs.cast(lsst.afw.image.makeWcs(md))
        result.append(tanwcs)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repo", metavar="repo", help="Directory containing butler repository.")
    parser.add_argument("refcat", metavar="refcat", help="Directory of reference catalog to match against.")
    parser.add_argument("ra", default=None, type=float, nargs='?',
                        help="Decimal degrees center RA of catalog (default uses center of first visit).")
    parser.add_argument("dec", default=None, type=float, nargs='?',
                        help="Decimal degrees center Dec of catalog (default uses center of first visit).")
    parser.add_argument("radius", default=6., type=float, nargs='?',
                        help="Radius (degrees) of sources to load from reference catalog.")
    parser.add_argument("-i", "--interactive", action="store_true",
                        help="Use interactive matplotlib backend and set ion(), in addition to saving files.")
    parser.add_argument("-o", "--outdir", default=".plots",
                        help="output directory for plots (default: $(default)s)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print extra things during calculations.")
    args = parser.parse_args()

    butler = lsst.daf.persistence.Butler(inputs=args.repo)
    dataIds = get_valid_dataIds(butler)

    data_refs = [butler.dataRef('wcs', dataId=dataId) for dataId in dataIds]
    visits = [data_ref.dataId['visit'] for data_ref in data_refs]
    old_wcs_list = get_old_wcs_list(data_refs)

    os.environ['ASTROMETRY_NET_DATA_DIR'] = args.refcat
    if args.ra is None or args.dec is None:
        calexp = data_refs[0].get('calexp')
        center = calexp.getInfo().getVisitInfo().getBoresightRaDec()
    else:
        center = (args.ra, args.dec)

    reference = prep_reference_loader(center, args.radius*degrees)

    jointcalStatistics = utils.JointcalStatistics(verbose=args.verbose)
    jointcalStatistics.compute_rms(data_refs, visits, reference)

    name = os.path.basename(os.path.normpath(args.repo))
    jointcalStatistics.make_plots(data_refs, visits, old_wcs_list, name=name,
                                  interactive=args.interactive, outdir=args.outdir)


if __name__ == "__main__":
    main()
