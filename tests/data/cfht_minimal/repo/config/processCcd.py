import lsst.pipe.tasks.processCcd
assert type(config)==lsst.pipe.tasks.processCcd.ProcessCcdConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.processCcd.ProcessCcdConfig' % (type(config).__module__, type(config).__name__)
import lsst.obs.cfht.cfhtIsrTask
# Save calibration results?
config.calibrate.doWrite=True

# Run deblender input exposure
config.calibrate.doDeblend=True

# Write reference matches (ignored if doWrite false)?
config.calibrate.doWriteMatches=True

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. Ignored if icSourceCat is not provided.
config.calibrate.icSourceFieldsToCopy=['calib_psfCandidate', 'calib_psfUsed', 'calib_psfReserved']

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.calibrate.requireAstrometry=True

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
# 	Valid Range = [0.0,inf)
config.calibrate.photoCal.magErrFloor=0.0

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.calibrate.photoCal.applyColorTerms=True

# use median instead of mean to compute zeropoint
config.calibrate.photoCal.useMedian=True

# Write a field name astrom_usedByPhotoCal to the schema
config.calibrate.photoCal.doWriteOutput=True

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.calibrate.photoCal.fluxField='slot_CalibFlux_flux'

config.calibrate.photoCal.colorterms.data={}
config.calibrate.photoCal.colorterms.data['e2v']=lsst.pipe.tasks.colorterms.ColortermDict()
config.calibrate.photoCal.colorterms.data['e2v'].data={}
config.calibrate.photoCal.colorterms.data['e2v'].data['u']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['u'].c2=0.0

# name of primary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['u'].primary='u'

# First-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['u'].c1=-0.241

# Constant parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['u'].c0=0.0

# name of secondary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['u'].secondary='g'

config.calibrate.photoCal.colorterms.data['e2v'].data['i']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['i'].c2=0.0

# name of primary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['i'].primary='i'

# First-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['i'].c1=0.085

# Constant parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['i'].c0=0.0

# name of secondary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['i'].secondary='r'

config.calibrate.photoCal.colorterms.data['e2v'].data['r']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['r'].c2=0.0

# name of primary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['r'].primary='r'

# First-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['r'].c1=0.024

# Constant parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['r'].c0=0.0

# name of secondary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['r'].secondary='g'

config.calibrate.photoCal.colorterms.data['e2v'].data['z']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['z'].c2=0.0

# name of primary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['z'].primary='z'

# First-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['z'].c1=-0.074

# Constant parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['z'].c0=0.0

# name of secondary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['z'].secondary='i'

config.calibrate.photoCal.colorterms.data['e2v'].data['g']=lsst.pipe.tasks.colorterms.Colorterm()
# Second-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['g'].c2=0.0

# name of primary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['g'].primary='g'

# First-order parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['g'].c1=-0.153

# Constant parameter
config.calibrate.photoCal.colorterms.data['e2v'].data['g'].c0=0.0

# name of secondary filter
config.calibrate.photoCal.colorterms.data['e2v'].data['g'].secondary='r'

# Use the extendedness parameter to select objects to use in photometric calibration?
# This applies only to the sources detected on the exposure, not the reference catalog
config.calibrate.photoCal.doSelectUnresolved=True

# number of iterations
config.calibrate.photoCal.nIter=20

# List of source flag fields that must be set for a source to be used.
config.calibrate.photoCal.goodFlags=[]

# Don't use objects fainter than this magnitude
config.calibrate.photoCal.magLimit=22.0

# clip at nSigma
config.calibrate.photoCal.nSigma=3.0

# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.calibrate.photoCal.photoCatName='e2v'

# maximum sigma to use when clipping
config.calibrate.photoCal.sigmaMax=0.25

# List of source flag fields that will cause a source to be rejected when they are set.
config.calibrate.photoCal.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# correction factor for modelFlux error
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for psfFlux error
config.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.calibrate.catalogCalculation.plugins.names=['base_ClassificationExtendedness']
# Scale factor to apply to shape for aperture
config.calibrate.measurement.undeblended['base_Variance'].scale=5.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Variance'].doMeasure=True

# Mask planes to ignore
config.calibrate.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Convergence tolerance for FWHM
config.calibrate.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# Additional value to add to background
config.calibrate.measurement.undeblended['base_SdssShape'].background=0.0

# Convergence tolerance for e1,e2
config.calibrate.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.calibrate.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum centroid shift, limited to 2-10
config.calibrate.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Maximum number of iterations
config.calibrate.measurement.undeblended['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SdssShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# fiddle factor for adjusting the binning
config.calibrate.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.calibrate.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.calibrate.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_GaussianCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_GaussianCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# Radius (in pixels) of apertures.
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_PixelFlags'].doMeasure=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.calibrate.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.calibrate.measurement.undeblended['base_Blendedness'].doShape=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.calibrate.measurement.undeblended['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.calibrate.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# FIXME! NEVER DOCUMENTED!
config.calibrate.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.measurement.undeblended['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.undeblended['base_FPPosition'].doMeasure=True

config.calibrate.measurement.undeblended.names=[]
# Scale factor to apply to shape for aperture
config.calibrate.measurement.plugins['base_Variance'].scale=5.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Mask planes to ignore
config.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Convergence tolerance for FWHM
config.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Additional value to add to background
config.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Convergence tolerance for e1,e2
config.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.calibrate.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum centroid shift, limited to 2-10
config.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Maximum number of iterations
config.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.5

# fiddle factor for adjusting the binning
config.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_GaussianCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# Radius (in pixels) of apertures.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.calibrate.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.calibrate.measurement.plugins['base_Blendedness'].doShape=True

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.calibrate.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.calibrate.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# FIXME! NEVER DOCUMENTED!
config.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.calibrate.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.calibrate.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

config.calibrate.measurement.plugins.names=['base_SdssShape', 'base_PsfFlux', 'base_SdssCentroid', 'base_GaussianCentroid', 'base_SkyCoord', 'base_CircularApertureFlux', 'base_PixelFlags', 'base_Variance', 'base_GaussianFlux', 'base_NaiveCentroid']
# Prefix to give undeblended plugins
config.calibrate.measurement.undeblendedPrefix='undeblended_'

# the name of the centroiding algorithm used to set source x,y
config.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source aperture flux slot
config.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the flux measurement algorithm used for calibration
config.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source psf flux slot
config.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source model flux slot
config.calibrate.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source inst flux slot
config.calibrate.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.calibrate.measurement.slots.shape='base_SdssShape'

# When measuring, replace other detected footprints with noise?
config.calibrate.measurement.doReplaceWithNoise=True

# Add ann offset to the generated noise.
config.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	variance	Mean = 0, variance = the image's variance
# 	measure	Measure clipped mean and variance from the whole image
# 
config.calibrate.measurement.noiseReplacer.noiseSource='measure'

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Run subtask to apply aperture correction
config.calibrate.doApCorr=True

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.calibrate.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.maxFootprintArea=1000000

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.calibrate.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.calibrate.deblend.psfChisq2b=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.calibrate.deblend.psfChisq1=1.5

# Guarantee that all peaks produce a child source.
config.calibrate.deblend.propagateAllPeaks=False

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	None	Field is optional
# 
config.calibrate.deblend.strayFluxRule='trim'

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.calibrate.deblend.weightTemplates=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.calibrate.deblend.catchFailures=False

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.minFootprintAxisRatio=0.0

# Find stray flux---flux not claimed by any child in the deblender.
config.calibrate.deblend.findStrayFlux=True

# Mask planes to ignore when performing statistics
config.calibrate.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.calibrate.deblend.maxTempDotProd=0.5

# When splitting stray flux, clip fractions below this value to zero.
config.calibrate.deblend.clipStrayFluxFraction=0.001

# Try to remove similar templates?
config.calibrate.deblend.removeDegenerateTemplates=False

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.calibrate.deblend.assignStrayFlux=True

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.calibrate.deblend.psfChisq2=1.5

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	ramp	Ramp down flux at the image edge by the PSF
# 	None	Field is optional
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 
config.calibrate.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	None	Field is optional
# 	necessary	When there is not an extended object in the footprint
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 
config.calibrate.deblend.strayFluxToPointSources='necessary'

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.calibrate.deblend.maskLimits={}

# Mask name for footprints not deblended, or None
config.calibrate.deblend.notDeblendedMask='NOT_DEBLENDED'

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.calibrate.requirePhotoCal=True

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.calibrate.astrometry.wcsFitter.numIter=3

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.order=3

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.calibrate.astrometry.wcsFitter.rejSigma=3.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.numRejIter=1

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.minMatchedPairs=30

# maximum determinant of linear transformation matrix for a usable solution
config.calibrate.astrometry.matcher.maxDeterminant=0.02

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceFluxType='Ap'

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.calibrate.astrometry.matcher.maxRotationDeg=1.0

# number of points to define a shape for matching
config.calibrate.astrometry.matcher.numPointsForShape=6

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matcher.maxMatchDistArcSec=5.0

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.calibrate.astrometry.matcher.minSnr=40.0

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.calibrate.astrometry.matcher.maxOffsetPix=300

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.calibrate.astrometry.matcher.allowedNonperpDeg=3.0

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.calibrate.astrometry.maxIter=3

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.calibrate.astrometry.forceKnownWcs=False

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.minMatchDistanceArcSec=0.001

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matchDistanceSigma=2.0

# Do temporary interpolated background subtraction before footprint detection?
config.calibrate.detection.doTempLocalBackground=False

# Estimate the background again after final source detection?
config.calibrate.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.calibrate.detection.nSigmaToGrow=2.4

# specifies the desired flavor of Threshold
# Allowed values:
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 
config.calibrate.detection.thresholdType='stdev'

# Pixels should be grown as isotropically as possible (slower)
config.calibrate.detection.isotropicGrow=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.calibrate.detection.returnOriginalFootprints=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.calibrate.detection.includeThresholdMultiplier=1.0

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.calibrate.detection.minPixels=1

# Fiddle factor to add to the background; debugging only
config.calibrate.detection.adjustBackground=0.0

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	both	detect both positive and negative sources
# 	negative	detect only negative sources
# 
config.calibrate.detection.thresholdPolarity='positive'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.calibrate.detection.background.algorithm='NATURAL_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detection.background.binSize=128

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.calibrate.detection.background.statisticsProperty='MEANCLIP'

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.calibrate.detection.background.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.background.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.background.weighting=True

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.background.useApprox=True

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.calibrate.detection.thresholdValue=5.0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detection.tempLocalBackground.binSize=64

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Names of mask planes to ignore while estimating the background
config.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.calibrate.detection.tempLocalBackground.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detection.tempLocalBackground.weighting=True

# Use Approximate (Chebyshev) to model background.
config.calibrate.detection.tempLocalBackground.useApprox=False

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.calibrate.applyApCorr.ignoreList=[]

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.calibrate.applyApCorr.proxies={}

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.calibrate.applyApCorr.doFlagApCorrFailures=True

# Perform astrometric calibration?
config.calibrate.doAstrometry=True

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.calibrate.checkUnitsParseStrict='raise'

# Perform phometric calibration?
config.calibrate.doPhotoCal=True

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.calibrate.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.calibrate.refObjLoader.filterMap={'i2': 'i'}

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.calibrate.refObjLoader.pixelMargin=50

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.calibrate.matchRadiusPix=3.0

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.calibrate.doWriteHeavyFootprintsInSources=True

import lsst.obs.cfht.cfhtIsrTask
config.isr.retarget(target=lsst.obs.cfht.cfhtIsrTask.CfhtIsrTask, ConfigClass=lsst.obs.cfht.cfhtIsrTask.CfhtIsrTaskConfig)

# The saturation level to use if no Detector is present in the Exposure (ignored if NaN)
config.isr.saturation=float('nan')

# FWHM of PSF (arcsec)
config.isr.fwhm=1.0

# Assemble amp-level exposures into a ccd-level exposure?
config.isr.doAssembleCcd=True

# fields to remove from the metadata of the assembled ccd.
config.isr.keysToRemoveFromAssembledCcd=[]

# Threshold used to stop iterating the brighter fatter correction.  It is the  absolute value of the difference between the current corrected image and the one from the previous iteration summed over all the pixels.
config.isr.brighterFatterThreshold=1000.0

# Maximum number of iterations for the brighter fatter correction
config.isr.brighterFatterMaxIter=10

# The method for fitting the overscan bias level.
# Allowed values:
# 	CUBIC_SPLINE	Fit cubic spline to the longest axis of the overscan region
# 	MEAN	Correct using the mean of the overscan region
# 	None	Field is optional
# 	LEG	Fit Legendre polynomial to the longest axis of the overscan region
# 	MEDIAN	Correct using the median of the overscan region
# 	AKIMA_SPLINE	Fit Akima spline to the longest axis of the overscan region
# 	POLY	Fit ordinary polynomial to the longest axis of the overscan region
# 	NATURAL_SPLINE	Fit natural spline to the longest axis of the overscan region
# 	CHEB	Fit Chebyshev polynomial to the longest axis of the overscan region
# 
config.isr.overscanFitType='MEDIAN'

# Apply fringe correction?
config.isr.doFringe=False

# Dataset type for input data; users will typically leave this alone, but camera-specific ISR tasks will override it
config.isr.datasetType='raw'

# renormalize to a gain of 1? (ignored if setGain false). Setting to True gives 1 ADU per electron. Setting to True is not recommended for mosaic cameras because it breaks normalization across the focal plane. However, if the CCDs are sufficiently flat then the resulting error may be acceptable.
config.isr.assembleCcd.doRenorm=False

# FITS headers to remove (in addition to DATASEC, BIASSEC, TRIMSEC and perhaps GAIN)
config.isr.assembleCcd.keysToRemove=[]

# set gain?
config.isr.assembleCcd.setGain=False

# trim out non-data regions?
config.isr.assembleCcd.doTrim=True

# The method for scaling the flat on the fly.
# Allowed values:
# 	MEAN	Scale by the inverse of the mean
# 	USER	Scale by flatUserScale
# 	MEDIAN	Scale by the inverse of the median
# 	None	Field is optional
# 
config.isr.flatScalingType='USER'

# The approximate flux of a zero-magnitude object in a one-second exposure
config.isr.fluxMag0T1=10000000000.0

# Rejection threshold (sigma) for collapsing overscan before fit
config.isr.overscanRej=3.0

# Order of polynomial or to fit if overscan fit type is a polynomial, or number of spline knots if overscan fit type is a spline.
config.isr.overscanOrder=1

# Number of fitting iterations
config.isr.fringe.iterations=20

# Sigma clip threshold
config.isr.fringe.clip=3.0

# Half-size of small (fringe) measurements (pixels)
config.isr.fringe.small=1

# Offset to the random number generator seed (full seed includes exposure ID)
config.isr.fringe.stats.rngSeedOffset=0

# Number of fitting iterations
config.isr.fringe.stats.iterations=3

# Ignore pixels with these masks
config.isr.fringe.stats.badMaskPlanes=['SAT']

# Statistic to use
config.isr.fringe.stats.stat=32

# Sigma clip threshold
config.isr.fringe.stats.clip=3.0

# Remove fringe pedestal?
config.isr.fringe.pedestal=True

# Half-size of large (background) measurements (pixels)
config.isr.fringe.large=50

# Only fringe-subtract these filters
config.isr.fringe.filters=['i', 'i2', 'z']

# Number of fringe measurements
config.isr.fringe.num=30000

# Apply flat field correction?
config.isr.doFlat=False

# Correct for nonlinearity of the detector's response?
config.isr.doLinearize=True

# Do fringe subtraction after flat-fielding?
config.isr.fringeAfterFlat=False

# update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip
config.isr.setGainAssembledCcd=True

# Number of pixels by which to grow the saturation footprints
config.isr.growSaturationFootprintSize=1

# Fallback default filter name for calibrations
config.isr.fallbackFilterName=None

# The gain to use if no Detector is present in the Exposure (ignored if NaN)
config.isr.gain=float('nan')

# Safety margin for CFHT sensors gain determination
config.isr.safe=0.95

# The read noise to use if no Detector is present in the Exposure
config.isr.readNoise=0.0

# Apply the brighter fatter correction
config.isr.doBrighterFatter=False

# If flatScalingType is 'USER' then scale flat by this amount; ignored otherwise
config.isr.flatUserScale=1.0

# Name of mask plane to use in saturation detection and interpolation
config.isr.saturatedMaskName='SAT'

# Apply dark frame correction?
config.isr.doDark=False

# Assemble amp-level calibration exposures into ccd-level exposure?
config.isr.doAssembleIsrExposures=True

# Apply bias frame correction?
config.isr.doBias=False

# Name of mask plane to use for suspect pixels
config.isr.suspectMaskName='SUSPECT'

# Persist postISRCCD?
config.isr.doWrite=False

# Should the gain be applied when applying the brighter fatter correction?
config.isr.brighterFatterApplyGain=True

# Width and height of PSF model, in pixels. Must be odd.
# 	Valid Range = [1,inf)
config.charImage.installSimplePsf.width=11

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.charImage.installSimplePsf.fwhm=3.5322300675464238

# Run deblender input exposure
config.charImage.doDeblend=False

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.charImage.checkUnitsParseStrict='raise'

# Measure PSF? If False then keep the existing PSF model (which must exist) and use that model for all operations.
config.charImage.doMeasurePsf=True

# number of times to look for contaminated pixels near known CR pixels
config.charImage.repair.cosmicray.niteration=3

# Don't interpolate over CR pixels
config.charImage.repair.cosmicray.keepCRs=False

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac2=0.4

# CRs must have > this many DN (== electrons/gain) in initial detection
config.charImage.repair.cosmicray.min_DN=150.0

# CRs must be > this many sky-sig above sky
config.charImage.repair.cosmicray.minSigma=6.0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.repair.cosmicray.background.binSize=100000

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# Names of mask planes to ignore while estimating the background
config.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.charImage.repair.cosmicray.background.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.weighting=True

# Use Approximate (Chebyshev) to model background.
config.charImage.repair.cosmicray.background.useApprox=False

# maximum number of contaminated pixels
config.charImage.repair.cosmicray.nCrPixelMax=100000

# Interpolate over defects? (ignored unless you provide a list of defects)
config.charImage.repair.doInterpolate=True

# Find and mask out cosmic rays?
config.charImage.repair.doCosmicRay=True

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.minSize=5

# Add a Gaussian to represent wings?
config.charImage.repair.interp.modelPsf.addWing=True

# Default FWHM of Gaussian model of core of star (pixels)
config.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.charImage.repair.interp.modelPsf.sizeFactor=3.0

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.charImage.repair.interp.modelPsf.size=None

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# Smoothly taper to the fallback value at the edge of the image?
config.charImage.repair.interp.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.charImage.repair.interp.negativeFallbackAllowed=True

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.charImage.repair.interp.fallbackUserValue=0.0

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.charImage.repair.interp.fallbackValueType='MEANCLIP'

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.charImage.astrometry.wcsFitter.numIter=3

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.maxScatterArcsec=10.0

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.order=4

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.charImage.astrometry.wcsFitter.rejSigma=3.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.numRejIter=1

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.charImage.astrometry.matcher.numBrightStars=50

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.charImage.astrometry.matcher.minMatchedPairs=30

# maximum determinant of linear transformation matrix for a usable solution
config.charImage.astrometry.matcher.maxDeterminant=0.02

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.charImage.astrometry.matcher.minFracMatchedPairs=0.3

# Type of source flux; typically one of Ap or Psf
config.charImage.astrometry.matcher.sourceFluxType='Ap'

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.charImage.astrometry.matcher.maxRotationDeg=1.0

# number of points to define a shape for matching
config.charImage.astrometry.matcher.numPointsForShape=6

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.charImage.astrometry.matcher.maxMatchDistArcSec=3.0

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.charImage.astrometry.matcher.minSnr=40.0

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.charImage.astrometry.matcher.maxOffsetPix=300

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.charImage.astrometry.matcher.allowedNonperpDeg=3.0

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.charImage.astrometry.maxIter=3

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.charImage.astrometry.forceKnownWcs=False

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.charImage.astrometry.minMatchDistanceArcSec=0.001

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.charImage.astrometry.matchDistanceSigma=2.0

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf is True then 2 should be plenty; otherwise more may be wanted.
# 	Valid Range = [1,inf)
config.charImage.psfIterations=2

# correction factor for modelFlux error
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for psfFlux error
config.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.charImage.catalogCalculation.plugins.names=['base_ClassificationExtendedness']
# Number of iterations for sigma clipping
config.charImage.measureApCorr.numIter=4

# Number of standard devisations to clip at
config.charImage.measureApCorr.numSigmaClip=3.0

# maximum Chebyshev function order in y
config.charImage.measureApCorr.fitConfig.orderY=2

# maximum Chebyshev function order in x
config.charImage.measureApCorr.fitConfig.orderX=2

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.charImage.measureApCorr.fitConfig.triangular=True

# Field name prefix for the flux other measurements should be aperture corrected to match
config.charImage.measureApCorr.refFluxName='slot_CalibFlux'

# Allow these measurement algorithms to fail without an exception
config.charImage.measureApCorr.allowFailure=[]

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measureApCorr.starSelector['secondMoment'].fluxLim=12500.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['secondMoment'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# size of the kernel to create
config.charImage.measureApCorr.starSelector['secondMoment'].kernelSize=21

# Multiplier of mean for maximum moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.charImage.measureApCorr.starSelector['secondMoment'].clumpNSigma=2.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measureApCorr.starSelector['secondMoment'].fluxMax=0.0

# Clipping threshold for moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentClip=5.0

# Number of bins in moment histogram
config.charImage.measureApCorr.starSelector['secondMoment'].histSize=64

# Multiplier of mean for minimum moments histogram range
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# Maximum moment to consider
config.charImage.measureApCorr.starSelector['secondMoment'].histMomentMax=100.0

# Name of a flag field that is True for stars that should be used.
config.charImage.measureApCorr.starSelector['flagged'].field='calib_psfUsed'

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['flagged'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# size of the kernel to create
config.charImage.measureApCorr.starSelector['flagged'].kernelSize=21

# maximum width to include in histogram
config.charImage.measureApCorr.starSelector['objectSize'].widthMax=10.0

# Keep objects within this many sigma of cluster 0's median
config.charImage.measureApCorr.starSelector['objectSize'].nSigmaClip=2.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['objectSize'].borderWidth=0

# size of the kernel to create
config.charImage.measureApCorr.starSelector['objectSize'].kernelSize=21

# Standard deviation of width allowed to be interpreted as good stars
config.charImage.measureApCorr.starSelector['objectSize'].widthStdAllowed=0.15

# minimum width to include in histogram
config.charImage.measureApCorr.starSelector['objectSize'].widthMin=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measureApCorr.starSelector['objectSize'].fluxMax=0.0

# Name of field in Source to use for flux measurement
config.charImage.measureApCorr.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measureApCorr.starSelector['objectSize'].fluxMin=12500.0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.charImage.measureApCorr.starSelector['catalog'].fluxMax=0.0

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.charImage.measureApCorr.starSelector['catalog'].fluxLim=0.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measureApCorr.starSelector['catalog'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measureApCorr.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# size of the kernel to create
config.charImage.measureApCorr.starSelector['catalog'].kernelSize=21

config.charImage.measureApCorr.starSelector.name='flagged'
# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.charImage.measureApCorr.minDegreesOfFreedom=1

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.charImage.doWriteExposure=True

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.charImage.background.algorithm='NATURAL_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.background.binSize=128

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.charImage.background.statisticsProperty='MEANCLIP'

# Names of mask planes to ignore while estimating the background
config.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.charImage.background.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.background.weighting=True

# Use Approximate (Chebyshev) to model background.
config.charImage.background.useApprox=True

# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.charImage.useSimplePsf=True

# Scale factor to apply to shape for aperture
config.charImage.measurement.undeblended['base_Variance'].scale=5.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Variance'].doMeasure=True

# Mask planes to ignore
config.charImage.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Convergence tolerance for FWHM
config.charImage.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# Additional value to add to background
config.charImage.measurement.undeblended['base_SdssShape'].background=0.0

# Convergence tolerance for e1,e2
config.charImage.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.charImage.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum centroid shift, limited to 2-10
config.charImage.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Maximum number of iterations
config.charImage.measurement.undeblended['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SdssShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.charImage.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.charImage.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# fiddle factor for adjusting the binning
config.charImage.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.charImage.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.charImage.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_GaussianCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_GaussianCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# Radius (in pixels) of apertures.
config.charImage.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.charImage.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_PixelFlags'].doMeasure=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.charImage.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.charImage.measurement.undeblended['base_Blendedness'].doShape=True

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.charImage.measurement.undeblended['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.charImage.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# Scaling factor of PSF FWHM for aperture radius.
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# FIXME! NEVER DOCUMENTED!
config.charImage.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# Value to subtract from the image pixel values
config.charImage.measurement.undeblended['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.charImage.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.undeblended['base_FPPosition'].doMeasure=True

config.charImage.measurement.undeblended.names=[]
# Scale factor to apply to shape for aperture
config.charImage.measurement.plugins['base_Variance'].scale=5.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Variance'].doMeasure=True

# Mask planes to ignore
config.charImage.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.charImage.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Convergence tolerance for FWHM
config.charImage.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Additional value to add to background
config.charImage.measurement.plugins['base_SdssShape'].background=0.0

# Convergence tolerance for e1,e2
config.charImage.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Whether to also compute the shape of the PSF model
config.charImage.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum centroid shift, limited to 2-10
config.charImage.measurement.plugins['base_SdssShape'].maxShift=0.0

# Maximum number of iterations
config.charImage.measurement.plugins['base_SdssShape'].maxIter=100

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SdssShape'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.charImage.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.charImage.measurement.plugins['base_Jacobian'].pixelScale=0.5

# fiddle factor for adjusting the binning
config.charImage.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.charImage.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SdssCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.charImage.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_GaussianCentroid'].doFootprintCheck=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_GaussianCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PeakCentroid'].doMeasure=True

# Radius (in pixels) of apertures.
config.charImage.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.charImage.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_PixelFlags'].doMeasure=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.charImage.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.charImage.measurement.plugins['base_Blendedness'].doShape=True

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.charImage.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.charImage.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# Scaling factor of PSF FWHM for aperture radius.
config.charImage.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# FIXME! NEVER DOCUMENTED!
config.charImage.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_GaussianFlux'].doMeasure=True

# Value to subtract from the image pixel values
config.charImage.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Do check that the centroid is contained in footprint.
config.charImage.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.charImage.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.charImage.measurement.plugins['base_FPPosition'].doMeasure=True

config.charImage.measurement.plugins.names=['base_SdssShape', 'base_PsfFlux', 'base_SdssCentroid', 'base_CircularApertureFlux', 'base_PixelFlags', 'base_GaussianFlux']
# Prefix to give undeblended plugins
config.charImage.measurement.undeblendedPrefix='undeblended_'

# the name of the centroiding algorithm used to set source x,y
config.charImage.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source aperture flux slot
config.charImage.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the flux measurement algorithm used for calibration
config.charImage.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source psf flux slot
config.charImage.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source model flux slot
config.charImage.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source inst flux slot
config.charImage.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.charImage.measurement.slots.shape='base_SdssShape'

# When measuring, replace other detected footprints with noise?
config.charImage.measurement.doReplaceWithNoise=True

# Add ann offset to the generated noise.
config.charImage.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	variance	Mean = 0, variance = the image's variance
# 	measure	Measure clipped mean and variance from the whole image
# 
config.charImage.measurement.noiseReplacer.noiseSource='measure'

# The seed multiplier value to use for random number generation
#    >= 1: set the seed deterministically based on exposureId
#       0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.charImage.measurement.noiseReplacer.noiseSeedMultiplier=1

# Do temporary interpolated background subtraction before footprint detection?
config.charImage.detection.doTempLocalBackground=False

# Estimate the background again after final source detection?
config.charImage.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.charImage.detection.nSigmaToGrow=2.4

# specifies the desired flavor of Threshold
# Allowed values:
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 
config.charImage.detection.thresholdType='stdev'

# Pixels should be grown as isotropically as possible (slower)
config.charImage.detection.isotropicGrow=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.charImage.detection.returnOriginalFootprints=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.charImage.detection.includeThresholdMultiplier=10.0

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.charImage.detection.minPixels=1

# Fiddle factor to add to the background; debugging only
config.charImage.detection.adjustBackground=0.0

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	both	detect both positive and negative sources
# 	negative	detect only negative sources
# 
config.charImage.detection.thresholdPolarity='positive'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.charImage.detection.background.algorithm='NATURAL_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.background.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detection.background.binSize=128

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.charImage.detection.background.statisticsProperty='MEANCLIP'

# Names of mask planes to ignore while estimating the background
config.charImage.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.charImage.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.charImage.detection.background.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.background.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detection.background.weighting=True

# Use Approximate (Chebyshev) to model background.
config.charImage.detection.background.useApprox=True

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.charImage.detection.thresholdValue=5.0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 	NONE	No background estimation is to be attempted
# 
config.charImage.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.approxOrderY=-1

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detection.tempLocalBackground.binSize=64

# type of statistic to use for grid points
# Allowed values:
# 	MEAN	unclipped mean
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	None	Field is optional
# 
config.charImage.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Names of mask planes to ignore while estimating the background
config.charImage.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 
config.charImage.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# Ignore NaNs when estimating the background
config.charImage.detection.tempLocalBackground.isNanSafe=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.approxOrderX=6

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detection.tempLocalBackground.weighting=True

# Use Approximate (Chebyshev) to model background.
config.charImage.detection.tempLocalBackground.useApprox=False

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measurePsf.starSelector['secondMoment'].fluxLim=12500.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['secondMoment'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['secondMoment'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter']

# size of the kernel to create
config.charImage.measurePsf.starSelector['secondMoment'].kernelSize=21

# Multiplier of mean for maximum moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMaxMultiplier=5.0

# candidate PSF's shapes must lie within this many sigma of the average shape
config.charImage.measurePsf.starSelector['secondMoment'].clumpNSigma=2.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measurePsf.starSelector['secondMoment'].fluxMax=0.0

# Clipping threshold for moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentClip=5.0

# Number of bins in moment histogram
config.charImage.measurePsf.starSelector['secondMoment'].histSize=64

# Multiplier of mean for minimum moments histogram range
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMinMultiplier=2.0

# Maximum moment to consider
config.charImage.measurePsf.starSelector['secondMoment'].histMomentMax=100.0

# Name of a flag field that is True for stars that should be used.
config.charImage.measurePsf.starSelector['flagged'].field='calib_psfUsed'

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['flagged'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['flagged'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# size of the kernel to create
config.charImage.measurePsf.starSelector['flagged'].kernelSize=21

# maximum width to include in histogram
config.charImage.measurePsf.starSelector['objectSize'].widthMax=10.0

# Keep objects within this many sigma of cluster 0's median
config.charImage.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['objectSize'].borderWidth=0

# size of the kernel to create
config.charImage.measurePsf.starSelector['objectSize'].kernelSize=21

# Standard deviation of width allowed to be interpreted as good stars
config.charImage.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# minimum width to include in histogram
config.charImage.measurePsf.starSelector['objectSize'].widthMin=0.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measurePsf.starSelector['objectSize'].fluxMax=0.0

# Name of field in Source to use for flux measurement
config.charImage.measurePsf.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_flux'

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measurePsf.starSelector['objectSize'].fluxMin=12500.0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
# 	Valid Range = [0.0,inf)
config.charImage.measurePsf.starSelector['catalog'].fluxMax=0.0

# specify the minimum psfFlux for good Psf Candidates
# 	Valid Range = [0.0,inf)
config.charImage.measurePsf.starSelector['catalog'].fluxLim=0.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector['catalog'].borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector['catalog'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter']

# size of the kernel to create
config.charImage.measurePsf.starSelector['catalog'].kernelSize=21

config.charImage.measurePsf.starSelector.name='objectSize'
# This number will be multiplied by the exposure ID to set the random seed for reserving candidates
config.charImage.measurePsf.reserveSeed=1

# Fraction of PSF candidates to reserve from fitting; none if <= 0
config.charImage.measurePsf.reserveFraction=-1.0

# Maximum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# size of cell used to determine PSF (pixels, column direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# specify spatial order for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# Reject candidates that are blended?
config.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# number of stars per psf cell for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# size of cell used to determine PSF (pixels, row direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# Minimum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# Rejection threshold (stdev) for candidates based on spatial fit
config.charImage.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# Mask blends in image?
config.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

# Use non-linear fitter for spatial variation of Kernel
config.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# for psf candidate evaluation
config.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# number of iterations of PSF candidate star list
config.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# Should each PSF candidate be given the same weight, independent of magnitude?
config.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# tolerance of spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# number of stars per psf Cell for spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# floor for variance is lam*data
config.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# number of eigen components for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

config.charImage.measurePsf.psfDeterminer.name='pca'
# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.charImage.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.maxFootprintArea=1000000

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.charImage.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.charImage.deblend.psfChisq2b=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.charImage.deblend.psfChisq1=1.5

# Guarantee that all peaks produce a child source.
config.charImage.deblend.propagateAllPeaks=False

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	None	Field is optional
# 
config.charImage.deblend.strayFluxRule='trim'

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.charImage.deblend.weightTemplates=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.charImage.deblend.catchFailures=False

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.minFootprintAxisRatio=0.0

# Find stray flux---flux not claimed by any child in the deblender.
config.charImage.deblend.findStrayFlux=True

# Mask planes to ignore when performing statistics
config.charImage.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.charImage.deblend.maxTempDotProd=0.5

# When splitting stray flux, clip fractions below this value to zero.
config.charImage.deblend.clipStrayFluxFraction=0.001

# Try to remove similar templates?
config.charImage.deblend.removeDegenerateTemplates=False

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.charImage.deblend.assignStrayFlux=True

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.charImage.deblend.psfChisq2=1.5

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	ramp	Ramp down flux at the image edge by the PSF
# 	None	Field is optional
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 
config.charImage.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	None	Field is optional
# 	necessary	When there is not an extended object in the footprint
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 
config.charImage.deblend.strayFluxToPointSources='necessary'

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.charImage.deblend.maskLimits={}

# Mask name for footprints not deblended, or None
config.charImage.deblend.notDeblendedMask='NOT_DEBLENDED'

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.charImage.applyApCorr.ignoreList=[]

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.charImage.applyApCorr.proxies={}

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.charImage.applyApCorr.doFlagApCorrFailures=True

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.charImage.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.charImage.refObjLoader.filterMap={}

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.charImage.refObjLoader.pixelMargin=50

# Run subtasks to measure and apply aperture corrections
config.charImage.doApCorr=True

# Persist results?
config.charImage.doWrite=True

# Perform calibration?
config.doCalibrate=True

