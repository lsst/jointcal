# This file is part of jointcal.

# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Functions to help create jointcal tests by generating fake data."""

__all__ = ['createFakeCatalog', 'createTwoFakeCcdImages', 'getMeasuredStarsFromCatalog']

import os
import numpy as np

import lsst.afw.geom
import lsst.afw.table
import lsst.daf.persistence
import lsst.obs.lsstSim
import lsst.pipe.base

import lsst.jointcal.star


def createTwoFakeCcdImages(num1=4, num2=4, seed=100, fakeCcdId=12):
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
        Sensor identifier to use for both CcdImages. The wcs, bbox, calib, etc.
        will still be drawn from the CFHT ccd=12 files, as that is the only
        testdata that is included in this simple test dataset.

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
    """
    np.random.seed(seed)

    visit1 = 849375
    visit2 = 850587
    instFluxKeyName = "SomeFlux"

    # Load or fake the necessary metadata for each CcdImage
    dataDir = lsst.utils.getPackageDir('jointcal')
    inputDir = os.path.join(dataDir, 'tests/data/cfht_minimal')
    butler = lsst.daf.persistence.Butler(inputDir)

    # so we can access parts of the camera later (e.g. focal plane)
    camera = butler.get('camera', visit=visit1)

    struct1 = createFakeCcdImage(butler, visit1, num1, instFluxKeyName,
                                 photoCalibMean=100.0, photoCalibErr=1.0, fakeCcdId=fakeCcdId)
    struct2 = createFakeCcdImage(butler, visit2, num2, instFluxKeyName,
                                 photoCalibMean=120.0, photoCalibErr=5.0, fakeCcdId=fakeCcdId)

    return lsst.pipe.base.Struct(camera=camera,
                                 catalogs=[struct1.catalog, struct2.catalog],
                                 ccdImageList=[struct1.ccdImage, struct2.ccdImage],
                                 bbox=struct1.bbox)


def createFakeCcdImage(butler, visit, num, instFluxKeyName,
                       photoCalibMean=100.0, photoCalibErr=1.0, fakeCcdId=12):
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
    instFluxKeyName : `str`
        Name of the instFluxKey to populate in the catalog.
    photoCalibMean : `float`, optional
        Value to set for calibrationMean in the created PhotoCalib.
    photoCalibErr : `float`, optional
        Value to set for calibrationErr in the created PhotoCalib.
    fakeCcdId : `int`, optional
        Use this as the ccdId in the returned CcdImage.

    Returns
    -------
    struct : `lsst.pipe.base.Struct`
       Result struct with components:

       - `catalog` : Catalogs containing fake sources
           (`lsst.afw.table.SourceCatalog`).
       - `ccdImage` : CcdImage containing the metadata and fake sources
           (`lsst.jointcal.CcdImage`).
       - `bbox` : Bounding Box of the image (`lsst.afw.geom.Box2I`).
    """
    ccdId = 12  # we only have data for ccd=12

    dataId = dict(visit=visit, ccd=ccdId)
    skyWcs = butler.get('calexp_wcs', dataId=dataId)
    visitInfo = butler.get('calexp_visitInfo', dataId=dataId)
    bbox = butler.get('calexp_bbox', dataId=dataId)
    detector = butler.get('calexp_detector', dataId=dataId)
    filt = butler.get("calexp_filter", dataId=dataId).getName()
    photoCalib = lsst.afw.image.PhotoCalib(photoCalibMean, photoCalibErr)

    catalog = createFakeCatalog(num, bbox, instFluxKeyName, skyWcs=skyWcs)
    ccdImage = lsst.jointcal.ccdImage.CcdImage(catalog, skyWcs, visitInfo, bbox, filt, photoCalib,
                                               detector, visit, fakeCcdId, instFluxKeyName)

    return lsst.pipe.base.Struct(catalog=catalog, ccdImage=ccdImage, bbox=bbox)


def createFakeCatalog(num, bbox, instFluxKeyName, skyWcs=None, refCat=False):
    """Return a fake minimally-useful catalog for jointcal.

    Parameters
    ----------
    num : `int`
        Number of sources to put in the catalogs. Should be
        a square, to have sqrt(num) centroids on a grid.
    bbox : `lsst.afw.geom.Box2I`
        Bounding Box of the detector to populate.
    instFluxKeyName : `str`
        Name of the instFluxKey to populate in the catalog.
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
    xErrKey = schema.addField("centroid_xSigma", type="F")
    yErrKey = schema.addField("centroid_ySigma", type="F")
    # shape
    shapeKey = lsst.afw.table.QuadrupoleKey.addFields(schema, "shape", "",
                                                      lsst.afw.table.CoordinateType.PIXEL)
    # Put the fake sources in the minimal catalog.
    schema.addField(instFluxKeyName+"_flux", type="D", doc="post-ISR instFlux")
    schema.addField(instFluxKeyName+"_fluxSigma", type="D", doc="post-ISR instFlux stddev")
    schema.addField(instFluxKeyName+"_calFlux", type="D", doc="maggies")
    schema.addField(instFluxKeyName+"_calFluxErr", type="D", doc="maggies stddev")
    schema.addField(instFluxKeyName+"_mag", type="D", doc="magnitude")
    schema.addField(instFluxKeyName+"_magErr", type="D", doc="magnitude stddev")
    return fillCatalog(schema, num, bbox,
                       centroidKey, xErrKey, yErrKey, shapeKey, instFluxKeyName,
                       skyWcs=skyWcs, refCat=refCat)


def fillCatalog(schema, num, bbox,
                centroidKey, xErrKey, yErrKey, shapeKey, instFluxKeyName,
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
    instFluxKeyName : `str`
        Name of instFlux field to populate (i.e. instFluxKeyName+'_flux')
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
    table.defineInstFlux(instFluxKeyName)
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
    catalog[instFluxKeyName + '_flux'] = instFlux
    catalog[instFluxKeyName + '_fluxSigma'] = instFluxErr

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
        star.setInstFlux(record.getInstFlux())
        star.setInstFluxErr(record.getInstFluxErr())
        # TODO: cleanup after DM-4044
        point = lsst.afw.geom.Point2D(star.x, star.y)
        pointFocal = pixToFocal.applyForward(point)
        star.setXFocal(pointFocal.getX())
        star.setYFocal(pointFocal.getY())
        stars.append(star)

    return stars
