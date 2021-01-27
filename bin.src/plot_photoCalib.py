#!/usr/bin/env python
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

"""
Plot a photoCalib over a full focal-plane for multiple visits.

Writes `photoCalib-tract_TRACT-visit_VISIT*.png` files to the current directory.

The persisted photocalib, after the individual ccd scalings are removed, should
be continuous across CCD boundaries, smooth, and relatively flat. Any major
variations in it suggest something unpleasant going on with the fit. Watch for
divergence along the bounaries: there may not have been enough stars to properly
constrain it.

The plots that divide out the original `singleFrameTask` calibs (either the
mean, or one per ccd) will show how the fitted calibration compares with the
original calib(s).
"""

import numpy as np
import lsst.afw.cameraGeom.utils
import lsst.daf.persistence

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.axes_grid1 import make_axes_locatable  # noqa: E402

matplotlib.rcParams['figure.dpi'] = 300


def getValidDataIds(butler, tract, dataset_type='jointcal_photoCalib'):
    """Return a list of all dataIds that exist in this butler.

    This exists here because the butler doesn't provide this functionality.

    Parameters
    ----------
    butler : `lsst.daf.persistence.Butler`
        The butler to search within.
    tract : `int`
        Tract id to include in search.
    dataset_type : `str`
        Dataset type to search for.
    """
    data_ids = []
    for data_ref in butler.subset(dataset_type, tract=tract):
        if data_ref.datasetExists():
            data_ids.append(data_ref.dataId)
    return data_ids


def colorbar(mappable):
    """Create a colorbar that obeys tight_layout, etc.

    Stolen from: http://joseph-long.com/writing/colorbars/

    Parameters
    ----------
    mappable
        Return value from e.g. `matplotlib.imshow()` or `matplotlib.scatter()`.
    """
    try:
        # QuadContourSet is not a Artist, so doesn't have `axes`
        ax = mappable.axes
    except AttributeError:
        ax = mappable.ax
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def makePhotoCalibImages(visit, butler, step=8, chips=[], tract=None,
                         chipScaling=False, meanCalib=False, singleCalib=False,
                         verbose=False):
    """Return a dict of the effect of each ccd's jointcal calibration.

    Parameters
    ----------
    visit : `int`
        The visit id to plot.
    butler : `lsst.daf.persistence.Butler`
        The butler instance to pull the data from.
    step : `int`
        Step size between samples.
        i.e. `np.linspace(0, bbox.getHeight(), bbox.getHeight()/binSize)`
    chips : `list` of `int`
        The chips to inc
    chipScaling : `bool`
        Include the constant per-chip component of the model.
        This will result in discontinuities across chip boundaries.
    meanCalib : `bool`
        Divide out the mean of the singleFrame `Calib` from each visit.
    singleCalib : `bool`
        Divide out each exposure's singleFrame `Calib`.
    verbose : `bool`
        Print statistics as we make the plots.
    """
    images = {}

    if meanCalib:
        meanCalibScaling = 0
        for ccd in chips:
            calib = butler.get('calexp_photoCalib', dataId=dict(visit=int(visit), ccd=int(ccd), tract=tract))
            meanCalibScaling += calib.getFluxMag0()[0]
        meanCalibScaling /= len(chips)
        if verbose:
            print('calib mean: %s' % meanCalibScaling)
    else:
        meanCalibScaling = 1

    for ccd in chips:
        detector = butler.get('calexp_detector', ccd=int(ccd))
        bbox = detector.getBBox()

        xx = np.linspace(bbox.getMinX(), bbox.getMaxX(), int(np.floor(bbox.getWidth()/step)))
        yy = np.linspace(bbox.getMinY(), bbox.getMaxY(), int(np.floor(bbox.getHeight()/step)))
        xgrid, ygrid = np.meshgrid(xx, yy)

        calibScaling = 1
        if singleCalib:
            calib = butler.get('calexp_photoCalib', dataId=dict(visit=int(visit), ccd=int(ccd), tract=tract))
            calibScaling = calib.getFluxMag0()[0]
            if verbose:
                print('calib ccd %s: %s' % (ccd, calibScaling))

        photoCalib = butler.get('jointcal_photoCalib', dataId=dict(visit=int(visit),
                                                                   ccd=int(ccd),
                                                                   tract=tract))
        if verbose:
            print("photoCalib mean ccd %s: %s" % (ccd, photoCalib.getCalibrationMean()))
        scaled = photoCalib.computeScaledCalibration()
        scaling = photoCalib.getCalibrationMean() if chipScaling else 1

        zgrid = scaled.evaluate(xgrid.flatten(), ygrid.flatten()).reshape(xgrid.shape)
        zgrid *= scaling * meanCalibScaling * calibScaling
        if verbose:
            print("mean of ccd %s: %s" % (ccd, zgrid.mean()))

        images[detector.getId()] = lsst.afw.image.Image(zgrid, dtype=zgrid.dtype.type)
    return images


class ImageMaker:
    """An image factory for lsst.afw.cameraGeom.utils.makeImageFromCamera().

    Inspired by lsst.pipe.drivers.makeCameraImage.
    """
    def __init__(self, images):
        self.isTrimmed = True
        self.images = images
        self.background = np.nan

    def getCcdImage(self, detector, imageFactory, binSize):
        """Return an image for a given detector."""
        if detector.getId() not in self.images:
            return None, detector
        else:
            return self.images[detector.getId()], detector


def plotVisitPhotoCalib(visit, butler, step=8, percentile=0,
                        tract=None, chips=[],
                        chipScaling=False, meanCalib=False, singleCalib=False,
                        colormap="magma", verbose=False):
    """Plot the post-jointcal calibration, and write it to a png.

    Parameters
    ----------
    visit : `int`
        The visit id to plot.
    butler : `lsst.daf.persistence.Butler`
        The butler instance to pull the data from.
    step : `int`
        Step size between samples.
        i.e. `np.linspace(0, bbox.getHeight(), bbox.getHeight()/binSize)`
    percentile : `float`
        Plot the data from percentile:100-percentile (chop off extreme outliers).
    chipScaling : `bool`
        Include the constant per-chip component of the model.
        This will result in discontinuities across chip boundaries.
    meanCalib : `bool`
        Divide out the mean of the singleFrame `Calib`s from the each exposure.
    singleCalib : `bool`
        Divide out each exposure's singleFrame `Calib` from each exposure.
    colormap : `str`
        matplotlib.Colormap name to use when making the plot.
    verbose : `bool`
        Print statistics as we make the plots.

    """
    filename = 'photoCalib-tract_%d-visit_%d'%(tract, visit)
    title = "visit: %d"%visit
    if chipScaling:
        filename += '-chipScale'
        title += " (chipScale)"
    if meanCalib:
        filename += '-meanCalib'
        title += " (meanCalib)"
    if singleCalib:
        filename += '-singleCalib'
        title += " (singleCalib)"
    filename += '.png'

    imageMaker = ImageMaker(makePhotoCalibImages(visit,
                                                 butler,
                                                 step=step,
                                                 chips=chips,
                                                 tract=tract,
                                                 chipScaling=chipScaling,
                                                 meanCalib=meanCalib,
                                                 singleCalib=singleCalib,
                                                 verbose=verbose))
    image = lsst.afw.cameraGeom.utils.makeImageFromCamera(butler.get('camera'),
                                                          imageSource=imageMaker,
                                                          detectorNameList=chips,
                                                          imageFactory=lsst.afw.image.ImageD,
                                                          binSize=step
                                                          )
    if verbose:
        print("Mean, Median of entire visit: %s, %s" %
              (np.nanmean(image.getArray()), np.nanmedian(image.getArray())))

    bbox = image.getBBox()
    extent = (0, bbox.getWidth()*step, 0, image.getHeight()*step)

    plt.figure()
    percentiles = np.nanpercentile(image.getArray(), [percentile, 100-percentile])
    ax = plt.imshow(image.getArray(), cmap=colormap, origin='lower',
                    extent=extent, vmin=percentiles[0], vmax=percentiles[1])
    plt.title(title)
    plt.xlabel("x focal plane pixels")
    plt.ylabel("y focal plane pixels")
    colorbar(ax)
    plt.tight_layout()

    plt.savefig(filename)
    print('Wrote', filename)


def main():
    """Commandline entry point for plot_photoCalib.py."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repo", metavar="repo",
                        help="Directory to butler repository containing `jointcal_photoCalib` datasets.")
    parser.add_argument("--step", default=8, type=int,
                        help="Steps between samples per image (i.e. larger number = coarser image).")
    parser.add_argument("--tract", default=0, type=int,
                        help="Tract to plot (until gen3 butler can manage this automatically).")
    parser.add_argument("--percentile", type=float, default=1,
                        help="Show only the percentile/100-percentile range of the photoCalib.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print extra information while producing plots.")
    parser.add_argument("--colormap", default="magma",
                        help="matplotlib.Colormap name to use when generating the 2d plot.")
    args = parser.parse_args()

    butler = lsst.daf.persistence.Butler(inputs=args.repo)
    dataIds = getValidDataIds(butler, args.tract)
    visits = np.unique([dataId['visit'] for dataId in dataIds])
    chips = np.unique([dataId['ccd'] for dataId in dataIds])
    if args.verbose:
        print("Found these visits:", visits)
        print("And these chips:", chips)

    for visit in visits:
        if args.verbose:
            print()
            print("Processing visit", visit)
        plotVisitPhotoCalib(visit, butler, tract=args.tract, chips=chips,
                            step=args.step, percentile=args.percentile,
                            chipScaling=False,
                            colormap=args.colormap, verbose=args.verbose)
        plotVisitPhotoCalib(visit, butler, tract=args.tract, chips=chips,
                            step=args.step, percentile=args.percentile,
                            chipScaling=True, meanCalib=True,
                            colormap=args.colormap, verbose=args.verbose)
        plotVisitPhotoCalib(visit, butler, tract=args.tract, chips=chips,
                            step=args.step, percentile=args.percentile,
                            chipScaling=True, singleCalib=True,
                            colormap=args.colormap, verbose=args.verbose)


if __name__ == "__main__":
    main()
