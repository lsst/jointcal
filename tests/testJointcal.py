# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os
import numpy as np
import collections

from lsst.afw import geom, coord, table
from lsst.meas.astrom import LoadAstrometryNetObjectsTask, LoadAstrometryNetObjectsConfig
import lsst.utils
from lsst.jointcal import jointcal

from lsst.daf import persistence


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
arcsec_per_radian = 206265.
absolute_error = 42e-3/arcsec_per_radian
# Set to True for a comparison plot and some diagnostic numbers.
doPlot = False


class FakeRef(object):
    """A mock data ref object, with minimal functionality."""

    def __init__(self, src, calexp_md, dataId):
        self.src = src
        self.calexp_md = calexp_md
        # info about our "telescope"
        # self.calexp_md.setDouble("AIRMASS", 1)
        self.calexp_md.setDouble("MJD", 54321.)
        # self.calexp_md.setDouble("EXPTIME", 30.)
        self.calexp_md.setDouble("LST", 53.00914)
        self.calexp_md.setDouble("RA2000", 53.00914)
        self.calexp_md.setDouble("DEC2000", -27.43895)

        self.dataId = dataId
        self.dataId['ccd'] = 1
        self.dataId['tract'] = 1

    def get(self, name, immediate=True):
        return getattr(self, name)

    def put(self, value, name):
        setattr(self, name, value)

    def getButler(self):
        class NamedThing(object):
            def __init__(self, name):
                self.name = name

            def getName(self):
                return self.name

        return {'camera': NamedThing('monkeySim')}


class JointcalTest(lsst.utils.tests.TestCase):
    def setUp(self):
        data_dir = lsst.utils.getPackageDir('validation_data_jointcal')
        os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'twinkles1_and_index')

        # position of the Twinkles run 1 catalog
        center = coord.IcrsCoord(53.00914*geom.degrees, -27.43895*geom.degrees)
        r = geom.Angle(3, geom.degrees)
        refLoader = LoadAstrometryNetObjectsTask(LoadAstrometryNetObjectsConfig())
        # Make a copy of the reference catalog for in-memory contiguity.
        self.reference = refLoader.loadSkyCircle(center, r, filterName='r').refCat.copy()

        # Get a set of source catalogs.
        self.catalogs = []
        butler = persistence.Butler(os.path.join(data_dir, 'twinkles1'))
        for visit in butler.queryMetadata('src', 'visit'):
            src = butler.get('src', dataId={'visit': visit})
            dataId = butler.dataRef('src', visit=visit).dataId
            calexp_md = butler.get('calexp_md', {'visit': visit, 'raft': '2,2', 'sensor': '1,1'})
            self.catalogs.append(FakeRef(src, calexp_md, dataId))

        self.jointcalTask = jointcal.JointcalTask()

    def _testJointCalTask(self, nCatalogs, relerr):
        """Test jointcal on nCatalogs, requiring < some relative error (arcsec)."""

        self.jointcalTask.run(self.catalogs[:nCatalogs])
        rms_rel, rms_abs = compute_rms(self.reference, self.catalogs[:nCatalogs], fluxlimit=100)
        # NOTE: this number come from an initial jointcal run.
        self.assertLess(rms_rel, relerr/arcsec_per_radian)
        self.assertLess(rms_abs, absolute_error)

    @unittest.skip('jointcal currently fails if only given one catalog!')
    def testJointCalTask_1_catalog(self):
        self.jointcalTask.run(self.catalogs[:1])

    def testJointCalTask_2_catalog(self):
        self._testJointCalTask(2, 8.4e-3)

    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_4_catalog(self):
        self._testJointCalTask(4, 7.8e-3)

    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_7_catalog(self):
        self._testJointCalTask(7, 7.5e-3)

    def testJointCalTask_10_catalog(self):
        self._testJointCalTask(10, 7.4e-3)


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


def make_match_dict(reference, catalogs, radius=0.1*geom.arcseconds, fluxlimit=100):
    """
    Return a dict of starID:[distances] over the catalogs, for RMS calculations.

    @param reference the reference catalog to use for absolute object positions
    @param catalogs (list) the catalogs to do the relative matching between
    @param radius (float arcseonds) source matching radius
    @param fluxlimit (float) SN (flux/fluxSigma) for sources to be included in match.

    @return (dict) dict of starID:list(measurement deltas for that star)
    """

    deltas = collections.defaultdict(list)
    for cat in catalogs:
        good = (cat.src.get('base_PsfFlux_flux')/cat.src.get('base_PsfFlux_fluxSigma')) > fluxlimit
        # things the classifier called stars.
        good &= (cat.src.get('base_ClassificationExtendedness_value') == 0)
        matches = table.matchRaDec(reference, cat.src[good], radius)
        for m in matches:
            deltas[m[0].getId()].append(m[2])
    return deltas


# NOTE: Stolen from meas_astrom.fitTanSipWcs.
def updateSourceCoords(wcs, catalog):
    """Update coords in a collection of sources, given a WCS."""
    schema = catalog.src.schema
    srcCoordKey = table.CoordKey(schema["coord"])
    for src in catalog.src:
        src.set(srcCoordKey, wcs.pixelToSky(src.getCentroid()))


def compute_rms(reference, catalogs, fluxlimit=100):
    """
    Match all catalogs to compute the RMS, for all detections above some flux limit.

    @param reference the reference catalog to use for absolute object positions
    @param catalogs (list) the catalogs to do the relative matching between
    @param fluxlimit (float) SN (flux/fluxSigma) for sources to be included in match.

    @return (float, float) relative, absolute RMS of stars, in radians
    """

    old_relative = make_match_dict(catalogs[0].src, catalogs[1:])
    old_absolute = make_match_dict(reference, catalogs, fluxlimit=fluxlimit)
    # Update coordinates with the new distortion corrections.
    for cat in catalogs:
        updateSourceCoords(cat.wcs.getWcs(), cat)
    new_relative = make_match_dict(catalogs[0].src, catalogs[1:])
    new_absolute = make_match_dict(reference, catalogs)

    if doPlot:
        old_rms_relative = rms_per_star(old_relative)
        old_rms_absolute = rms_per_star(old_absolute)
        new_rms_relative = rms_per_star(new_relative)
        new_rms_absolute = rms_per_star(new_absolute)
        print(len(catalogs))
        print("relative RMS:", rms_total(old_relative)*arcsec_per_radian,
              rms_total(new_relative)*arcsec_per_radian)
        print("absolute RMS:", rms_total(old_absolute)*arcsec_per_radian,
              rms_total(new_absolute)*arcsec_per_radian)
        plot_deltas(old_rms_relative, old_rms_absolute, new_rms_relative, new_rms_absolute)

    return rms_total(new_relative), rms_total(new_absolute)


def plot_deltas(old_rms_relative, old_rms_absolute, new_rms_relative, new_rms_absolute):
    import matplotlib.pyplot as plt
    plt.ion()

    plt.figure()

    color_rel = 'green'
    color_abs = 'blue'
    plt.title('relative vs. absolute: %d vs. %d'%(len(old_rms_relative), len(old_rms_absolute)))

    plt.hist(old_rms_relative*arcsec_per_radian, color=color_rel, ls='dashed', label='old rel',
             histtype='step', bins=30, lw=1.5)
    plt.hist(new_rms_relative*arcsec_per_radian, color=color_rel, ls='solid', label='new rel',
             histtype='step', bins=30, lw=1.5)

    plt.hist(old_rms_absolute*arcsec_per_radian, color=color_abs, ls='dashed', label='old abs',
             histtype='step', bins=30, lw=1.5)
    plt.hist(new_rms_absolute*arcsec_per_radian, color=color_abs, ls='solid', label='new abs',
             histtype='step', bins=30, lw=1.5)

    plt.xlabel('arcseconds')
    plt.legend(loc='best')

    # So one can muck-about with things after plotting...
    import ipdb
    ipdb.set_trace()


# TODO: the memory test cases currently fail in jointcal. I'll have to clean that up later.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
