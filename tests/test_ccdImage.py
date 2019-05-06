# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Test creation and use of the CcdImage class."""
import unittest

import lsst.utils.tests
from lsst.jointcal import testUtils

import lsst.afw.geom
import lsst.jointcal


class CcdImageTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.nStars1 = 4
        self.nStars2 = 100
        struct = testUtils.createTwoFakeCcdImages(num1=self.nStars1, num2=self.nStars2)
        self.ccdImage1 = struct.ccdImageList[0]
        self.ccdImage2 = struct.ccdImageList[1]
        self.bbox = struct.bbox

        self.associations = lsst.jointcal.Associations(struct.ccdImageList)

    def checkCountStars(self, ccdImage, nStars):
        """Check that ccdImage.countStars() is correct under various conditions.

        Parameters
        ----------
        ccdImage : `lsst.jointcal.CcdImage`
            The ccdImage to test.
        nStars : `int`
            Number of stars in ccdImage's catalog.

        Notes
        -----
        Does not test the case where some ``measuredStars`` are not ``valid``,
        as there is no interface for modifying that the python level. To test
        that would require creating a `Fitter` and a fake outlier list and
        calling ``removeMeasOutliers`` and/or ``removeRefOutliers``, which
        cannot be easily made part of the python API.
        """
        # By default there are no stars because catalogForFit is uninitialized.
        measStars, refStars = ccdImage.countStars()
        self.assertEqual(measStars, 0)
        self.assertEqual(refStars, 0)

        # With no associations, the catalog should have exactly as many valid
        # measuredStars as were created.
        ccdImage.resetCatalogForFit()
        measStars, refStars = ccdImage.countStars()
        self.assertEqual(measStars, nStars)
        self.assertEqual(refStars, 0)

        # Cross match catalogs: there will still be no refcat matches.
        matchCut = 3.0 * lsst.afw.geom.arcseconds
        # There should be no fittedStars until we associate the catalogs.
        self.assertEqual(self.associations.fittedStarListSize(), 0)
        self.associations.computeCommonTangentPoint()
        self.associations.associateCatalogs(matchCut)
        # Confirm that every measuredStar (in both ccdImages) got a fittedStar associated to it.
        self.assertEqual(self.associations.fittedStarListSize(), self.nStars1 + self.nStars2)
        # measuredStars and refStars should be unchanged after association.
        self.assertEqual(measStars, nStars)
        self.assertEqual(refStars, 0)

        # Make a fake reference catalog; will match the catalog one-to-one.
        skyWcs = ccdImage.getReadWcs().getSkyWcs()
        self.refCat = testUtils.createFakeCatalog(nStars, self.bbox, "refFlux", skyWcs=skyWcs, refCat=True)
        # associate the reference stars
        self.associations.collectRefStars(self.refCat, matchCut, 'refFlux_instFlux', 0.1)
        measStars, refStars = ccdImage.countStars()
        self.assertEqual(measStars, nStars)
        self.assertEqual(refStars, nStars)

    def testCcdImage1(self):
        self.checkCountStars(self.ccdImage1, self.nStars1)

    def testCcdImage2(self):
        self.checkCountStars(self.ccdImage2, self.nStars2)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
