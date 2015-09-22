from __future__ import division, absolute_import
#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import numpy

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.coadd.utils as coaddUtils
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord

from lsst.afw.fits import FitsError
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask, SelectStruct
from lsst.coadd.utils import CoaddDataIdContainer
from lsst.pipe.tasks.getRepositoryData import DataRefListRunner
from lsst.meas.astrom.loadAstrometryNetObjects import LoadAstrometryNetObjectsTask
from lsst.meas.astrom import AstrometryNetDataConfig

from lsst.meas.simastrom.simastromLib import SimAstromControl, simAstrom, Associations, ProjectionHandler , AstromFit, SimplePolyModel, OneTPPerShoot

__all__ = ["SimAstromConfig", "SimAstromTask"]

class SimAstromConfig(pexConfig.Config):
    """Config for SimAstromTask
    """

# Keep this config parameter as a place holder
    doWrite = pexConfig.Field(
        doc = "persist SimAstrom output...",
        dtype = bool,
        default = True,
    )
    
    sourceFluxField = pexConfig.Field(
        doc = "Type of source flux",
        dtype = str,
        default = "base_CircularApertureFlux_5",   # base_CircularApertureFlux_17_0 in recent stack version 
    )

class SimAstromTask(pipeBase.CmdLineTask):
 
    ConfigClass = SimAstromConfig
    RunnerClass = DataRefListRunner
    _DefaultName = "simAstrom"
    
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
#        self.makeSubtask("select")

# We don't need to persist config and metadata at this stage. In this way, we don't need to put a specific entry in the
# camera mapper policy file
    def _getConfigName(self):
        return None
        
    def _getMetadataName(self):
        return None
        
    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)

        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9")
        return parser

    @pipeBase.timeMethod
    def run(self, ref):
        
#        sourceCat = test2()
        
        
        configSel = StarSelectorConfig()
        ss = StarSelector(configSel, self.config.sourceFluxField)
        
        print self.config.sourceFluxField
        astromControl = SimAstromControl()
        astromControl.sourceFluxField = self.config.sourceFluxField

        assoc = Associations()
        
        for dataRef in ref :
            src = dataRef.get("src", immediate=True)
            calexp = dataRef.get("calexp", immediate=True)
            tanwcs = afwImage.TanWcs.cast(calexp.getWcs())
            bbox = calexp.getBBox()
            md = dataRef.get("calexp_md", immediate=True)
            calib = afwImage.Calib(md)
            filt = calexp.getFilter().getName()
            
            newSrc = ss.select(src, calib)
            if len(newSrc) == 0 :
                print "no source selected in ", dataRef.dataId["visit"], dataRef.dataId["ccd"]
                continue
            print "%d sources selected in visit %d - ccd %d"%(len(newSrc), dataRef.dataId["visit"], dataRef.dataId["ccd"])
            
        # Should call a source selector here in order to send a list
        # of reasonable star to the fitter.
            assoc.AddImage(newSrc, tanwcs, md, bbox, filt, calib,
                           dataRef.dataId['visit'], dataRef.dataId['ccd'],
                           dataRef.getButler().mapper.getCameraName(), 
                           astromControl)
        
        matchCut = 3.0
        assoc.AssociateCatalogs(matchCut)
        
        # Use external reference catalogs handled by LSST stack mechanism
        # Get the bounding box overlapping all associated images
        bbox = assoc.GetRaDecBBox()
        center = afwCoord.Coord(bbox.getCenter(), afwGeom.degrees)
        corner = afwCoord.Coord(bbox.getMax(), afwGeom.degrees)
        radius = center.angularSeparation(corner).asRadians()
        
        # Get astrometry_net_data path
        anDir = lsst.utils.getPackageDir('astrometry_net_data')
        if anDir is None:
            raise RuntimeError("astrometry_net_data is not setup")

        andConfig = AstrometryNetDataConfig()
        andConfigPath = os.path.join(anDir, "andConfig.py")
        if not os.path.exists(andConfigPath):
            raise RuntimeError("astrometry_net_data config file \"%s\" required but not found" %andConfigPath)
        andConfig.load(andConfigPath)
        
        task = LoadAstrometryNetObjectsTask.ConfigClass()
        loader = LoadAstrometryNetObjectsTask(task)
        
        # Determine default filter associated to the catalog
        filt, mfilt = andConfig.magColumnMap.items()[0]
        print "Using", filt, "band for reference flux"

        refCat = loader.loadSkyCircle(center, afwGeom.Angle(radius, afwGeom.radians), filt).refCat
        print refCat.getSchema().getOrderedNames()
        
#        assoc.CollectRefStars(False) # To use USNO-A catalog 

        assoc.CollectLSSTRefStars(refCat)
        assoc.SelectFittedStars()
        assoc.DeprojectFittedStars() # required for AstromFit
        sky2TP = OneTPPerShoot(assoc.TheCcdImageList())
        spm = SimplePolyModel(assoc.TheCcdImageList(), sky2TP, True, 0)

        fit = AstromFit(assoc,spm)
        fit.Minimize("Distortions")
        fit.Minimize("Positions")
        fit.Minimize("Distortions Positions")

        for i in range(100): 
            nout = fit.RemoveOutliers(5.) # 5 sigma
            fit.Minimize("Distortions Positions")
            chi2 = fit.ComputeChi2()
            
            print chi2
            if (nout == 0) : break
        fit.MakeResTuple("res.list")



class StarSelectorConfig(pexConfig.Config):
    
    badFlags = pexConfig.ListField(
        doc = "List of flags which cause a source to be rejected as bad",
        dtype = str,
        default = [ "base_PixelFlags_flag_saturated", 
                    "base_PixelFlags_flag_cr",
                    "base_PixelFlags_flag_interpolated",
                    "base_SdssCentroid_flag",
                    "base_SdssShape_flag"],
    )

class StarSelector(object) :
    
    ConfigClass = StarSelectorConfig

    def __init__(self, config, sourceFluxField):
        """Construct a star selector
        
        @param[in] config: An instance of StarSelectorConfig
        """
        self.config = config
        self.sourceFluxField = sourceFluxField
    
    def select(self, srcCat, calib):
# Return a catalog containing only reasonnable stars

        schema = srcCat.getSchema()
        newCat = afwTable.SourceCatalog(schema)
        fluxKey = schema[self.sourceFluxField+"_flux"].asKey()
        fluxErrKey = schema[self.sourceFluxField+"_fluxSigma"].asKey()
        parentKey = schema["parent"].asKey()
        flagKeys = []
        for f in self.config.badFlags :
            key = schema[f].asKey()
            flagKeys.append(key)
        fluxFlagKey = schema[self.sourceFluxField+"_flag"].asKey()
        flagKeys.append(fluxFlagKey)
        
        for src in srcCat :
            # Reject galaxies
#            if src.get("base_ClassificationExtendedness_value") > 0.5 :
#                continue
            # Do not consider sources with bad flags
            for f in flagKeys :
                rej = 0
                if src.get(f) :
                    rej = 1
                    break
            if rej == 1 :
                continue
            # Reject negative flux
            flux = src.get(fluxKey)
            if flux < 0 :
                continue
            # Reject objects with too large magnitude
            fluxErr = src.get(fluxErrKey)
            mag, magErr = calib.getMagnitude(flux, fluxErr)
            if mag > 22.5 or magErr > 0.1 or flux/fluxErr < 10 :
                continue
            # Reject blends
            if src.get(parentKey) != 0 :
                continue
            footprint = src.getFootprint()
            if footprint is not None and len(footprint.getPeaks()) > 1 :
                continue
                
            newCat.append(src)
            
#        print len(srcCat), len(newCat)
        
        return newCat
