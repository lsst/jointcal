# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os

from lsst.afw import geom, coord
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('validation_data_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'cfht_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
absolute_error = 48e-3/jointcalTestBase.arcsec_per_radian
# Set to True for a comparison plot and some diagnostic numbers.
doPlot = False


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHT(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        jointcalTestBase.JointcalTestBase.setUp(self)
        self.doPlot = doPlot
        self.matchRadius = 0.1*geom.arcseconds

        # position of the validation_data_cfht catalog
        center = coord.IcrsCoord(214.8848320609*geom.degrees, 52.6622198737*geom.degrees)
        radius = geom.Angle(3, geom.degrees)
        self._prep_reference_loader(center, radius)

        self.input_dir = os.path.join(data_dir, 'cfht')
        self.visitList = [849375, 850587]

    @unittest.skipIf(data_dir is None, "validation_data_jointcal not setup")
    def test_jointcalTask_2_visits(self):
        self._testJointCalTask(2, 25e-3, absolute_error)


# TODO: the memory test cases currently fail in jointcal. I'll have to clean that up later.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
