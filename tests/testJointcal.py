# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os
import numpy as np

import eups

from lsst.daf.base import PropertyList
from lsst.afw import geom, coord, table
from lsst.meas.astrom import astrometry
import lsst.utils
from lsst.jointcal import jointcal


class FakeRef(object):
    """A mock data ref object, with minimal functionality."""
    def __init__(self, refCat):
        """Generate a fake data ref, perturbed from a reference catalog."""

        # make the fake catalog, adding necessary fields.
        cat = refCat.copy(deep=True)
        # TODO: don't think I need to pre-select the stars here: jointcal is supposed to do that.
        # cat = cat[cat['resolved'] == False]  # noqa E712 is incorrect here!
        # randomly perturb the points.
        error = geom.Angle(1e-2, geom.arcseconds)
        # NOTE: have to use asRadians() since cat['coord*'] is a numpy view and doesn't know units.
        cat['coord_ra'] += np.random.normal(size=len(cat)) * error.asRadians() / np.cos(cat['coord_dec'])
        cat['coord_dec'] += np.random.normal(size=len(cat)) * error.asRadians()

        # make a new catalog, adding the new fields we need for a source catalog + errors
        mapper = table.SchemaMapper(cat.schema)
        mapper.addMinimalSchema(table.SourceTable.makeMinimalSchema(), True)
        for x in cat.schema:
            mapper.addMapping(x.key, True)
        for f in ['u', 'g', 'r', 'i', 'z', 'y']:
            mapper.editOutputSchema().addField('{}_fluxSigma'.format(f), type='F', doc='{} error'.format(f))
            mapper.editOutputSchema().addField('{}_flag'.format(f), type='F', doc='{} flags'.format(f))
        for f in ["base_PixelFlags_flag_saturated",
                  "base_PixelFlags_flag_cr",
                  "base_PixelFlags_flag_interpolated",
                  "base_SdssCentroid_flag",
                  "base_SdssShape_flag"]:
            mapper.editOutputSchema().addField(f, type='I')

        newCat = table.SourceCatalog(mapper.getOutputSchema())
        newCat.extend(cat, mapper=mapper)
        self.src = newCat

        # make some metadata for a fake WCS.
        self.calexp_md = PropertyList()
        self.calexp_md.set("SIMPLE", "T")
        self.calexp_md.set("BITPIX", -32)
        self.calexp_md.set("NAXIS", 2)
        self.calexp_md.set("NAXIS1", 1024)
        self.calexp_md.set("NAXIS2", 1153)
        self.calexp_md.set("RADECSYS", 'ICRS')
        self.calexp_md.set("EQUINOX", 2000.)
        self.calexp_md.setDouble("CRVAL1", 53.12499999999999)
        self.calexp_md.setDouble("CRVAL2", -28.1)
        self.calexp_md.setDouble("CRPIX1", 1109.99981456774)
        self.calexp_md.setDouble("CRPIX2", 560.018167811613)
        self.calexp_md.set("CTYPE1", 'RA---SIN')
        self.calexp_md.set("CTYPE2", 'DEC--SIN')
        self.calexp_md.setDouble("CD1_1", 5.10808596133527E-05)
        self.calexp_md.setDouble("CD1_2", 1.85579539217196E-07)
        self.calexp_md.setDouble("CD2_2", -5.10281493481982E-05)
        self.calexp_md.setDouble("CD2_1", -8.27440751733828E-07)

        # make a fake dataId
        self.dataId = {}
        self.dataId['filter'] = 'r'
        self.dataId['visit'] = 1
        self.dataId['ccd'] = 1
        self.dataId['tract'] = 1

    def get(self, name, immediate=True):
        return getattr(self, name)

    def getButler(self):
        class NamedThing(object):
            def getName(self):
                return 'monkey'
        return {'camera': NamedThing()}


class JointcalTest(lsst.utils.tests.TestCase):
    @classmethod
    def setUpClass(cls):
        """Can only do the eups setup once, so have to do it here only."""
        # NOTE: this may be cheating? How else should I do this?
        eups.setup('astrometry_net_data', productRoot=os.path.join(os.environ['JOINTCAL_DIR'], 'tests/data/'))

    def setUp(self):
        # Generate a few random realizations of a catalog.
        center = coord.IcrsCoord('03:32:30', '-28:06:00')  # position of the Twinkles run 1 catalog
        r = geom.Angle(1, geom.degrees)
        astrometryTask = astrometry.AstrometryTask()
        catalog = astrometryTask.refObjLoader.loadSkyCircle(center, r, filterName='r')
        self.catalogs = []
        for _ in range(3):
            self.catalogs.append(FakeRef(catalog.refCat))

        # tweak the jointcal runner to the config we want.
        self.jointcalTask = jointcal.JointcalTask()
        self.jointcalTask.config.sourceFluxField = 'r'

    def testJointCalTask(self):
        self.jointcalTask.run(self.catalogs)


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
