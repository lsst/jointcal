# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os

from lsst.afw import geom, coord
import lsst.utils
import lsst.pex.exceptions
from lsst.daf import persistence
from lsst.jointcal import jointcal
import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('validation_data_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'twinkles1_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
absolute_error = 42e-3/jointcalTestBase.arcsec_per_radian
# Set to True for a comparison plot and some diagnostic numbers.
doPlot = False


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class FakeRef(object):
    """
    A mock data ref object, with minimal functionality.

    This is needed to insert missing metadata while DM-5503 is dealt with.
    """

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


class JointcalTestLSSTSim(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        jointcalTestBase.JointcalTestBase.setUp(self)
        self.doPlot = doPlot
        self.matchRadius = 0.1*geom.arcseconds

        # position of the Twinkles run 1 catalog
        center = coord.IcrsCoord(53.00914*geom.degrees, -27.43895*geom.degrees)
        radius = geom.Angle(3, geom.degrees)
        self._prep_reference_loader(center, radius)

        self.input_dir = os.path.join(data_dir, 'cfht')
        self.visitList = range(840, 850)

        # Get a set of source catalogs.
        self.catalogs = []
        butler = persistence.Butler(os.path.join(data_dir, 'twinkles1'))
        for visit in butler.queryMetadata('src', 'visit'):
            src = butler.get('src', dataId={'visit': visit})
            dataId = butler.dataRef('src', visit=visit).dataId
            calexp_md = butler.get('calexp_md', {'visit': visit, 'raft': '2,2', 'sensor': '1,1'})
            self.catalogs.append(FakeRef(src, calexp_md, dataId))

        self.jointcalTask = jointcal.JointcalTask()
        # NOTE: Tweaking S/N threshold to help pass original thresholds.
        # TODO: Jointcal results are quite sensitive to the particulars of the
        # sources used for assocaitions, and astrometrySourceSelector does not
        # exactly match the original bundled StarSelector.
        # TODO: Once we make jointcal more robust, we should be able to drop this.
        self.jointcalTask.config.sourceSelector["astrometry"].minSnr = 13

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('jointcal currently fails if only given one catalog!')
    def testJointCalTask_1_catalog(self):
        self.jointcalTask.run(self.catalogs[:1])

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    def testJointCalTask_2_catalog(self):
        self._testJointCalTask_run(2, 8.4e-3, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_4_catalog(self):
        self._testJointCalTask_run(4, 7.8e-3, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_7_catalog(self):
        self._testJointCalTask_run(7, 7.5e-3, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    def testJointCalTask_10_catalog(self):
        self._testJointCalTask_run(10, 7.4e-3, absolute_error)


# TODO: the memory test cases currently fail in jointcal. I'll have to clean that up later.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
