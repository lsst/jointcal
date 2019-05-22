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

"""Functions to help create jointcal tests by generating fake data."""

__all__ = ['createFakeCatalog', 'createTwoFakeCcdImages', 'getMeasuredStarsFromCatalog']

import os
import unittest

import numpy as np

import lsst.afw.geom
import lsst.afw.table
import lsst.daf.persistence
import lsst.pipe.base

import lsst.jointcal.star


def canRunTests():
    """Returns True if the necessary packages and files are available.

    We need ``obs_cfht`` to load the test/data/cfht_minimal dataset, which
    includes the metadata that is used to build the fake catalogs.
    """
    try:
        import lsst.obs.cfht  # noqa: F401
        return True
    except ImportError:
        return False


def createTwoFakeCcdImages(num1=4, num2=4, seed=100, fakeCcdId=12,
                           photoCalibMean1=1e-2, photoCalibMean2=1.2e-2,
                           fakeWcses=(None, None),
                           fakeVisitInfos=(None, None)):
    """Return two fake ccdImages built on CFHT Megacam metadata.

    If ``num1 == num2``, the catalogs will align on-sky so each source will
    have a match in the other catalog.

    This uses the butler dataset stored in `tests/data/cfht_minimal` to
    bootstrap the metadata.

    Parameters
    ----------
    num1, num2 : `int`, optional
        Number of sources to put in the first and second catalogs. Should be
        a square, to have sqrt(num) centroids on a grid.
    seed : `int`, optional
        Seed value for np.random.
    fakeCcdId : `int`, optional
        Sensor identifier to use for both CcdImages. The wcs, bbox, photoCalib, etc.
        will still be drawn from the CFHT ccd=12 files, as that is the only
        testdata that is included in this simple test dataset.
    photoCalibMean1, photoCalibMean2: `float`, optional
        The mean photometric calibration to pass to each ccdImage construction.
        Note: this value is 1/instFluxMag0, so it should be less than 1.
    fakeWcses : `list` [`lsst.afw.geom.SkyWcs`], optional
        The SkyWcses to use instead of the ones read from disk.
    fakeWcses : `list` [`lsst.afw.image.VisitInfo`], optional
        The VisitInfos to use instead of the ones read from disk.

    Returns
    -------
    struct : `lsst.pipe.base.Struct`
       Result struct with components:

       - `camera` : Camera representing these catalogs
           (`lsst.afw.cameraGeom.Camera`).
       - `catalogs` : Catalogs containing fake sources
           (`list` of `lsst.afw.table.SourceCatalog`).
       - `ccdImageList` : CcdImages containing the metadata and fake sources
           (`list` of `lsst.jointcal.CcdImage`).
       - `bbox` : Bounding Box of the image (`lsst.afw.geom.Box2I`).
       - 'fluxFieldName' : name of the instFlux field in the catalogs ('str').
    """
    if not canRunTests():
        msg = "Necessary packages not available to run tests that use the cfht_minimal dataset."
        raise unittest.SkipTest(msg)

    np.random.seed(seed)

    visit1 = 849375
    visit2 = 850587
    fluxFieldName = "SomeFlux"

    # Load or fake the necessary metadata for each CcdImage
    dataDir = lsst.utils.getPackageDir('jointcal')
    inputDir = os.path.join(dataDir, 'tests/data/cfht_minimal')
    butler = lsst.daf.persistence.Butler(inputDir)

    # so we can access parts of the camera later (e.g. focal plane)
    camera = butler.get('camera', visit=visit1)

    struct1 = createFakeCcdImage(butler, visit1, num1, fluxFieldName,
                                 photoCalibMean=photoCalibMean1, photoCalibErr=1.0, fakeCcdId=fakeCcdId,
                                 fakeWcs=fakeWcses[0], fakeVisitInfo=fakeVisitInfos[0])
    struct2 = createFakeCcdImage(butler, visit2, num2, fluxFieldName,
                                 photoCalibMean=photoCalibMean2, photoCalibErr=5.0, fakeCcdId=fakeCcdId,
                                 fakeWcs=fakeWcses[1], fakeVisitInfo=fakeVisitInfos[1])

    return lsst.pipe.base.Struct(camera=camera,
                                 catalogs=[struct1.catalog, struct2.catalog],
                                 ccdImageList=[struct1.ccdImage, struct2.ccdImage],
                                 bbox=struct1.bbox,
                                 skyWcs=[struct1.skyWcs, struct2.skyWcs],
                                 fluxFieldName=fluxFieldName)


def createFakeCcdImage(butler, visit, num, fluxFieldName,
                       photoCalibMean=1e-2, photoCalibErr=1.0, fakeCcdId=12,
                       fakeWcs=None, fakeVisitInfo=None):
    """Create a fake CcdImage by making a fake catalog.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        Butler to load metadata from.
    visit : `int`
        Visit identifier to build a butler dataId.
    num : `int`
        Number of sources to put in the catalogs. Should be
        a square, to have sqrt(num) centroids on a grid.
    fluxFieldName : `str`
        Name of the flux field to populate in the catalog, without `_instFlux`
        (e.g. "slot_CalibFlux").
    photoCalibMean : `float`, optional
        Value to set for calibrationMean in the created PhotoCalib.
        Note: this value is 1/instFluxMag0, so it should be less than 1.
    photoCalibErr : `float`, optional
        Value to set for calibrationErr in the created PhotoCalib.
    fakeCcdId : `int`, optional
        Use this as the ccdId in the returned CcdImage.
    fakeWcs : `lsst.afw.geom.SkyWcs`, optional
        A SkyWcs to use instead of one read from disk.
    fakeVisitInfo : `lsst.afw.image.VisitInfo`, optional
        A VisitInfo to use instead of one read from disk.

    Returns
    -------
    struct : `lsst.pipe.base.Struct`
       Result struct with components:

       - `catalog` : Catalogs containing fake sources
           (`lsst.afw.table.SourceCatalog`).
       - `ccdImage` : CcdImage containing the metadata and fake sources
           (`lsst.jointcal.CcdImage`).
       - `bbox` : Bounding Box of the image (`lsst.afw.geom.Box2I`).
       - `skyWcs` : SkyWcs of the image (`lsst.afw.geom.SkyWcs`).
    """
    ccdId = 12  # we only have data for ccd=12

    dataId = dict(visit=visit, ccd=ccdId)
    skyWcs = fakeWcs if fakeWcs is not None else butler.get('calexp_wcs', dataId=dataId)
    visitInfo = fakeVisitInfo if fakeVisitInfo is not None else butler.get('calexp_visitInfo', dataId=dataId)
    bbox = butler.get('calexp_bbox', dataId=dataId)
    detector = butler.get('calexp_detector', dataId=dataId)
    filt = butler.get("calexp_filter", dataId=dataId).getName()
    photoCalib = lsst.afw.image.PhotoCalib(photoCalibMean, photoCalibErr)

    catalog = createFakeCatalog(num, bbox, fluxFieldName, skyWcs=skyWcs)
    ccdImage = lsst.jointcal.ccdImage.CcdImage(catalog, skyWcs, visitInfo, bbox, filt, photoCalib,
                                               detector, visit, fakeCcdId, fluxFieldName)

    return lsst.pipe.base.Struct(catalog=catalog, ccdImage=ccdImage, bbox=bbox, skyWcs=skyWcs)


def createFakeCatalog(num, bbox, fluxFieldName, skyWcs=None, refCat=False):
    """Return a fake minimally-useful catalog for jointcal.

    Parameters
    ----------
    num : `int`
        Number of sources to put in the catalogs. Should be
        a square, to have sqrt(num) centroids on a grid.
    bbox : `lsst.afw.geom.Box2I`
        Bounding Box of the detector to populate.
    fluxFieldName : `str`
        Name of the flux field to populate in the catalog, without `_instFlux`
        (e.g. "slot_CalibFlux").
    skyWcs : `lsst.afw.geom.SkyWcs` or None, optional
        If supplied, use this to fill in coordinates from centroids.
    refCat : `bool`, optional
        Return a ``SimpleCatalog`` so that it behaves like a reference catalog?

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        A populated source catalog.
    """
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    # centroid
    centroidKey = lsst.afw.table.Point2DKey.addFields(schema, "centroid", "centroid", "pixels")
    xErrKey = schema.addField("centroid_xErr", type="F")
    yErrKey = schema.addField("centroid_yErr", type="F")
    # shape
    shapeKey = lsst.afw.table.QuadrupoleKey.addFields(schema, "shape", "",
                                                      lsst.afw.table.CoordinateType.PIXEL)
    # Put the fake sources in the minimal catalog.
    schema.addField(fluxFieldName+"_instFlux", type="D", doc="post-ISR instFlux")
    schema.addField(fluxFieldName+"_instFluxErr", type="D", doc="post-ISR instFlux stddev")
    schema.addField(fluxFieldName+"_flux", type="D", doc="source flux (nJy)")
    schema.addField(fluxFieldName+"_fluxErr", type="D", doc="flux stddev (nJy)")
    schema.addField(fluxFieldName+"_mag", type="D", doc="magnitude")
    schema.addField(fluxFieldName+"_magErr", type="D", doc="magnitude stddev")
    return fillCatalog(schema, num, bbox,
                       centroidKey, xErrKey, yErrKey, shapeKey, fluxFieldName,
                       skyWcs=skyWcs, refCat=refCat)


def fillCatalog(schema, num, bbox,
                centroidKey, xErrKey, yErrKey, shapeKey, fluxFieldName,
                skyWcs=None, fluxErrFraction=0.05, refCat=False):
    """Return a catalog populated with fake, but reasonable, sources.

    Centroids are placed on a uniform grid, errors are normally distributed.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Pre-built schema to make the catalog from.
    num : `int`
        Number of sources to put in the catalog.
    bbox : `lsst.afw.geom.Box2I`
        Bounding box of the ccd to put sources in.
    centroidKey : `lsst.afw.table.Key`
        Key for the centroid field to populate.
    xErrKey : `lsst.afw.table.Key`
        Key for the xErr field to populate.
    yErrKey : `lsst.afw.table.Key`
        Key for the yErr field to populate.
    shapeKey : `lsst.afw.table.Key`
        Key for the shape field to populate.
    fluxFieldName : `str`
        Name of the flux field to populate in the catalog, without `_instFlux`
        (e.g. "slot_CalibFlux").
    skyWcs : `lsst.afw.geom.SkyWcs` or None, optional
        If supplied, use this to fill in coordinates from centroids.
    fluxErrFraction : `float`, optional
        Fraction of instFlux to use for the instFluxErr.
    refCat : `bool`, optional
        Return a ``SimpleCatalog`` so that it behaves like a reference catalog?

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`
        The filled catalog.
    """
    table = lsst.afw.table.SourceTable.make(schema)
    table.defineCentroid('centroid')
    table.defineShape('shape')
    table.defineCalibFlux(fluxFieldName)
    if refCat:
        catalog = lsst.afw.table.SimpleCatalog(table)
    else:
        catalog = lsst.afw.table.SourceCatalog(table)

    instFlux = np.random.random(num)*10000
    instFluxErr = np.abs(instFlux * np.random.normal(fluxErrFraction, scale=0.1, size=num))
    xx = np.linspace(bbox.getMinX(), bbox.getMaxX(), int(np.sqrt(num)))
    yy = np.linspace(bbox.getMinY(), bbox.getMaxY(), int(np.sqrt(num)))
    xv, yv = np.meshgrid(xx, yy)
    vx = np.random.normal(scale=0.1, size=num)
    vy = np.random.normal(scale=0.1, size=num)

    # make all the sources perfectly spherical, for simplicity.
    mxx = 1
    myy = 1
    mxy = 0

    for i, (x, y) in enumerate(zip(xv.ravel(), yv.ravel())):
        record = catalog.addNew()
        record.set('id', i)
        record.set(centroidKey, lsst.afw.geom.Point2D(x, y))
        record.set(shapeKey, lsst.afw.geom.ellipses.Quadrupole(mxx, myy, mxy))

    if skyWcs is not None:
        lsst.afw.table.updateSourceCoords(skyWcs, catalog)

    catalog[xErrKey] = vx
    catalog[yErrKey] = vy
    catalog[fluxFieldName + '_instFlux'] = instFlux
    catalog[fluxFieldName + '_instFluxErr'] = instFluxErr

    return catalog


def getMeasuredStarsFromCatalog(catalog, pixToFocal):
    """Return a list of measuredStars built from a catalog.

    Parameters
    ----------
    catalog : `lsst.afw.table.SourceCatalog`
        The table to get sources from.
    pixToFocal : `lsst.afw.geom.TransformPoint2ToPoint2`
        Transform that goes from pixel to focal plane coordinates, to set the
        MeasuredStar x/y focal points.

    Returns
    -------
    stars : `list` of `lsst.jointcal.MeasuredStar`
        MeasuredStars built from the catalog sources.
    """
    stars = []
    for record in catalog:
        star = lsst.jointcal.star.MeasuredStar()
        star.x = record.getX()
        star.y = record.getY()
        star.setInstFluxAndErr(record.getCalibInstFlux(), record.getCalibInstFluxErr())
        # TODO: cleanup after DM-4044
        point = lsst.afw.geom.Point2D(star.x, star.y)
        pointFocal = pixToFocal.applyForward(point)
        star.setXFocal(pointFocal.getX())
        star.setYFocal(pointFocal.getY())
        stars.append(star)

    return stars
