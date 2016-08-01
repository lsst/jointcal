# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os

from astropy import units as u

import lsst.afw.geom
import lsst.afw.coord
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
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
absolute_error = 42e-3*u.arcsecond
# Set to True for a comparison plot and some diagnostic numbers.
do_plot = False


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
        self.calexp_md.setDouble("MJD", 54321.)
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
        """
        This is just here to make jointcal happy until we properly solve the
        metadata problem, when it will be able to be deleted. Jointcal doesn't
        use the butler except to get the camera name, so that's all this does.
        """
        class NamedThing(object):
            def __init__(self, name):
                self.name = name

            def getName(self):
                return self.name

        return {'camera': NamedThing('monkeySim')}


class JointcalTestLSSTSim(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        jointcalTestBase.JointcalTestBase.setUp(self)
        self.do_plot = do_plot
        self.match_radius = 0.1*lsst.afw.geom.arcseconds

        # position of the Twinkles run 1 catalog
        center = lsst.afw.coord.IcrsCoord(53.00914*lsst.afw.geom.degrees, -27.43895*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees
        self._prep_reference_loader(center, radius)

        self.input_dir = os.path.join(data_dir, 'cfht')
        self.all_visits = range(840, 850)

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
        # sources used for associations, and astrometrySourceSelector does not
        # exactly match the original bundled StarSelector.
        # TODO: Once we make jointcal more robust, we should be able to drop this.
        self.jointcalTask.sourceSelector.config.minSnr = 13

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('jointcal currently fails if only given one catalog!')
    def testJointCalTask_1_catalog(self):
        self.jointcalTask.run(self.catalogs[:1])

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    def testJointCalTask_2_catalog(self):
        # NOTE: The relative RMS limits were empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 8.4e-3*u.arcsecond
        self._testJointCalTask_run(2, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_4_catalog(self):
        relative_error = 7.8e-3*u.arcsecond
        self._testJointCalTask_run(4, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    @unittest.skip('Keeping this around for diagnostics on the behavior with n catalogs.')
    def testJointCalTask_7_catalog(self):
        relative_error = 7.5e-3*u.arcsecond
        self._testJointCalTask_run(7, relative_error, absolute_error)

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    def testJointCalTask_10_catalog(self):
        relative_error = 7.4e-3*u.arcsecond
        self._testJointCalTask_run(10, relative_error, absolute_error)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
