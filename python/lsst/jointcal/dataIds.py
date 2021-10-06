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

# Copied from HyperSuprime-Cam/pipe_tasks
import collections

import lsst.log
import lsst.pex.exceptions
import lsst.afw.table
import lsst.afw.image
import lsst.meas.base

from lsst.coadd.utils import CoaddDataIdContainer

__all__ = ["PerTractCcdDataIdContainer"]


class PerTractCcdDataIdContainer(CoaddDataIdContainer):
    """A version of lsst.pipe.base.DataIdContainer that combines raw data IDs (defined as whatever we use
    for 'src') with a tract.
    """

    def castDataIds(self, butler):
        """Validate data IDs and cast them to the correct type (modify idList in place).

        Parameters
        ----------
        butler
            data butler
        """
        try:
            idKeyTypeDict = butler.getKeys(datasetType="src", level=self.level)
        except KeyError as e:
            raise KeyError("Cannot get keys for datasetType %s at level %s: %s" % ("src", self.level, e))

        idKeyTypeDict = idKeyTypeDict.copy()
        idKeyTypeDict["tract"] = int

        for dataDict in self.idList:
            for key, strVal in dataDict.items():
                try:
                    keyType = idKeyTypeDict[key]
                except KeyError:
                    validKeys = sorted(idKeyTypeDict.keys())
                    raise KeyError("Unrecognized ID key %r; valid keys are: %s" % (key, validKeys))
                if keyType != str:
                    try:
                        castVal = keyType(strVal)
                    except Exception:
                        raise TypeError("Cannot cast value %r to %s for ID key %r" % (strVal, keyType, key,))
                    dataDict[key] = castVal

    def _addDataRef(self, namespace, dataId, tract):
        """Construct a dataRef based on dataId, but with an added tract key"""
        forcedDataId = dataId.copy()
        forcedDataId['tract'] = tract
        dataRef = namespace.butler.dataRef(datasetType=self.datasetType, dataId=forcedDataId)
        self.refList.append(dataRef)

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        if self.datasetType is None:
            raise RuntimeError("Must call setDatasetType first")
        skymap = None
        log = lsst.log.Log.getLogger("jointcal.dataIds")
        visitTract = {}  # Set of tracts for each visit
        visitRefs = {}  # List of data references for each visit
        for dataId in self.idList:
            if "tract" not in dataId:
                # Discover which tracts the data overlaps
                log.infof("Reading WCS to determine tracts for components of dataId={}", dict(dataId))
                if skymap is None:
                    skymap = self.getSkymap(namespace)

                for ref in namespace.butler.subset("calexp", dataId=dataId):
                    if not ref.datasetExists("calexp"):
                        log.warning("calexp with dataId: %s not found.", dict(dataId))
                        continue

                    # XXX fancier mechanism to select an individual exposure than just pulling out "visit"?
                    if "visit" in ref.dataId.keys():
                        visit = ref.dataId["visit"]
                    else:
                        # Fallback if visit is not in the dataId
                        visit = namespace.butler.queryMetadata("calexp", ("visit"), ref.dataId)[0]
                    if visit not in visitRefs:
                        visitRefs[visit] = list()
                    visitRefs[visit].append(ref)

                    wcs = ref.get("calexp_wcs", immediate=True)
                    box = lsst.geom.Box2D(ref.get("calexp_bbox"))
                    # Going with just the nearest tract.  Since we're throwing all tracts for the visit
                    # together, this shouldn't be a problem unless the tracts are much smaller than a CCD.
                    tract = skymap.findTract(wcs.pixelToSky(box.getCenter()))
                    if lsst.meas.base.imageOverlapsTract(tract, wcs, box):
                        if visit not in visitTract:
                            visitTract[visit] = set()
                        visitTract[visit].add(tract.getId())
            else:
                tract = dataId.pop("tract")
                # making a DataRef for src fills out any missing keys and allows us to iterate
                for ref in namespace.butler.subset("src", dataId=dataId):
                    if ref.datasetExists():
                        self._addDataRef(namespace, ref.dataId, tract)

        # Ensure all components of a visit are kept together by putting them all in the same set of tracts
        # NOTE: sorted() here is to keep py2 and py3 dataRefs in the same order.
        # NOTE: see DM-9393 for details.
        for visit, tractSet in sorted(visitTract.items()):
            for ref in visitRefs[visit]:
                for tract in sorted(tractSet):
                    self._addDataRef(namespace, ref.dataId, tract)
        if visitTract:
            tractCounter = collections.Counter()
            for tractSet in visitTract.values():
                tractCounter.update(tractSet)
            log.infof("Number of visits per tract: {}", dict(tractCounter))
