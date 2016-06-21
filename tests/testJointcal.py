# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os
import numpy as np
import collections

import eups

from lsst.afw import geom, coord, table
from lsst.meas import astrom
import lsst.utils
from lsst.jointcal import jointcal

from lsst.daf import persistence

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
absolute_error = 42e-3/206265
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
    @classmethod
    def setUpClass(cls):
        """Can only do the eups setup once, so have to do it here only."""

        # NOTE: this may be cheating? How else should I do this?
        jointcal_reference_dir = os.path.join(os.environ['VALIDATION_DATA_JOINTCAL_DIR'],
                                              'twinkles1_and_index')
        eups.setup('astrometry_net_data', productRoot=jointcal_reference_dir)

    def setUp(self):

        data_dir = os.environ['VALIDATION_DATA_JOINTCAL_DIR']

        # position of the Twinkles run 1 catalog
        center = coord.IcrsCoord(53.00914*geom.degrees, -27.43895*geom.degrees)
        r = geom.Angle(3, geom.degrees)
        config = astrom.LoadAstrometryNetObjectsTask.ConfigClass()
        refLoader = astrom.LoadAstrometryNetObjectsTask(config=config)
        reference = refLoader.loadSkyCircle(center, r, filterName='r')
        self.reference = reference.refCat.copy()  # for in-memory contiguity.

        # Get a set of source catalogs.
        self.catalogs = []
        butler = persistence.Butler(os.path.join(data_dir, 'twinkles1'))
        for visit in butler.queryMetadata('src', 'visit'):
            src = butler.get('src', dataId={'visit': visit})
            dataId = butler.dataRef('src', visit=visit).dataId
            calexp_md = butler.get('calexp_md', {'visit': visit, 'raft': '2,2', 'sensor': '1,1'})
            self.catalogs.append(FakeRef(src, calexp_md, dataId))

        self.jointcalTask = jointcal.JointcalTask()

    @unittest.skip('jointcal currently fails if only given one catalog!')
    def testJointCalTask_1_catalog(self):
        self.jointcalTask.run(self.catalogs[:1])

    def testJointCalTask_2_catalog(self):
        ncat = 2
        self.jointcalTask.run(self.catalogs[:ncat])
        rms_rel, rms_abs = compute_rms(self.reference, self.catalogs[:ncat], fluxlimit=100)
        # NOTE: this number come from an initial jointcal run.
        self.assertLess(rms_rel, 8.4e-3/206265)
        self.assertLess(rms_abs, absolute_error)

    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_4_catalog(self):
        ncat = 4
        self.jointcalTask.run(self.catalogs[:ncat])
        rms_rel, rms_abs = compute_rms(self.reference, self.catalogs[:ncat], fluxlimit=100)
        # NOTE: this number come from an initial jointcal run.
        self.assertLess(rms_rel, 7.8e-3/206265)
        self.assertLess(rms_abs, absolute_error)

    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_7_catalog(self):
        ncat = 7
        self.jointcalTask.run(self.catalogs[:ncat])
        rms_rel, rms_abs = compute_rms(self.reference, self.catalogs[:ncat], fluxlimit=100)
        # NOTE: this number come from an initial jointcal run.
        self.assertLess(rms_rel, 7.5e-3/206265)
        self.assertLess(rms_abs, absolute_error)

    def testJointCalTask_10_catalog(self):
        ncat = 10
        self.jointcalTask.run(self.catalogs[:ncat])
        rms_rel, rms_abs = compute_rms(self.reference, self.catalogs[:ncat], fluxlimit=100)
        # NOTE: this number comes from an initial jointcal run.
        self.assertLess(rms_rel, 7.4e-3/206265)
        self.assertLess(rms_abs, absolute_error)


def rms_per_star(data):
    """Compute the RMS/catalog star; data is a dict of lists of distances."""
    result = np.empty(len(data))
    for i, k in enumerate(data):
        result[i] = np.sqrt(np.mean(np.array(data[k])**2))
    return result


def rms_total(data):
    """Compute the total rms across all stars."""
    n = 0
    total = 0
    for i, k in enumerate(data):
        total += sum(np.array(data[k])**2)
        n += len(data[k])
    return np.sqrt(total/n)


def make_match_dict(reference, catalogs, radius=0.1*geom.arcseconds, fluxlimit=100):
    """Make a dict of object:delta of all the cross-catalog matches, to calculate the RMS across."""

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
    """Match all catalogs to compute the RMS, for all detections above some flux limit."""

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
        print("relative RMS:", rms_total(old_relative)*206265, rms_total(new_relative)*206265)
        print("absolute RMS:", rms_total(old_absolute)*206265, rms_total(new_absolute)*206265)
        plot_deltas(old_rms_relative, old_rms_absolute, new_rms_relative, new_rms_absolute)

    return rms_total(new_relative), rms_total(new_absolute)


def plot_deltas(old_rms_relative, old_rms_absolute, new_rms_relative, new_rms_absolute):
    import matplotlib.pyplot as plt
    plt.ion()

    plt.figure()

    color_rel = 'green'
    color_abs = 'blue'
    plt.title('relative vs. absolute: %d vs. %d'%(len(old_rms_relative), len(old_rms_absolute)))

    # 206265 arcseconds per radian
    plt.hist(old_rms_relative*206265, color=color_rel, ls='dashed', label='old rel',
             histtype='step', bins=30, lw=1.5)
    plt.hist(new_rms_relative*206265, color=color_rel, ls='solid', label='new rel',
             histtype='step', bins=30, lw=1.5)

    plt.hist(old_rms_absolute*206265, color=color_abs, ls='dashed', label='old abs',
             histtype='step', bins=30, lw=1.5)
    plt.hist(new_rms_absolute*206265, color=color_abs, ls='solid', label='new abs',
             histtype='step', bins=30, lw=1.5)

    plt.xlabel('arcseconds')
    plt.legend(loc='best')

    # So one can muck-about with things after plotting...
    import ipdb
    ipdb.set_trace()


def suite():
    lsst.utils.tests.init()
    suites = []
    suites += unittest.makeSuite(JointcalTest)
    # suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
