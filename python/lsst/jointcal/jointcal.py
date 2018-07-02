# See COPYRIGHT file at the top of the source tree.
import collections
import numpy as np

import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.pex.exceptions as pexExceptions
import lsst.afw.table
import lsst.meas.algorithms
from lsst.verify import Job, Measurement

from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .dataIds import PerTractCcdDataIdContainer

import lsst.jointcal
from lsst.jointcal import MinimizeResult

__all__ = ["JointcalConfig", "JointcalRunner", "JointcalTask"]

Photometry = collections.namedtuple('Photometry', ('fit', 'model'))
Astrometry = collections.namedtuple('Astrometry', ('fit', 'model', 'sky_to_tan_projection'))


# TODO: move this to MeasurementSet in lsst.verify per DM-12655.
def add_measurement(job, name, value):
    meas = Measurement(job.metrics[name], value)
    job.measurements.insert(meas)


class JointcalRunner(pipeBase.ButlerInitializedTaskRunner):
    """Subclass of TaskRunner for jointcalTask

    jointcalTask.run() takes a number of arguments, one of which is a list of dataRefs
    extracted from the command line (whereas most CmdLineTasks' run methods take
    single dataRef, are are called repeatedly). This class transforms the processed
    arguments generated by the ArgumentParser into the arguments expected by
    Jointcal.run().

    See pipeBase.TaskRunner for more information.
    """

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """
        Return a list of tuples per tract, each containing (dataRefs, kwargs).

        Jointcal operates on lists of dataRefs simultaneously.
        """
        kwargs['profile_jointcal'] = parsedCmd.profile_jointcal
        kwargs['butler'] = parsedCmd.butler

        # organize data IDs by tract
        refListDict = {}
        for ref in parsedCmd.id.refList:
            refListDict.setdefault(ref.dataId["tract"], []).append(ref)
        # we call run() once with each tract
        result = [(refListDict[tract], kwargs) for tract in sorted(refListDict.keys())]
        return result

    def __call__(self, args):
        """
        Parameters
        ----------
        args
            Arguments for Task.run()

        Returns
        -------
        pipe.base.Struct
            if self.doReturnResults is False:

            - ``exitStatus``: 0 if the task completed successfully, 1 otherwise.

            if self.doReturnResults is True:

            - ``result``: the result of calling jointcal.run()
            - ``exitStatus``: 0 if the task completed successfully, 1 otherwise.
        """
        exitStatus = 0  # exit status for shell

        # NOTE: cannot call self.makeTask because that assumes args[0] is a single dataRef.
        dataRefList, kwargs = args
        butler = kwargs.pop('butler')
        task = self.TaskClass(config=self.config, log=self.log, butler=butler)
        result = None
        try:
            result = task.run(dataRefList, **kwargs)
            exitStatus = result.exitStatus
            job_path = butler.get('verify_job_filename')
            result.job.write(job_path[0])
        except Exception as e:  # catch everything, sort it out later.
            if self.doRaise:
                raise e
            else:
                exitStatus = 1
                eName = type(e).__name__
                tract = dataRefList[0].dataId['tract']
                task.log.fatal("Failed processing tract %s, %s: %s", tract, eName, e)

        if self.doReturnResults:
            return pipeBase.Struct(result=result, exitStatus=exitStatus)
        else:
            return pipeBase.Struct(exitStatus=exitStatus)


class JointcalConfig(pexConfig.Config):
    """Config for JointcalTask"""

    doAstrometry = pexConfig.Field(
        doc="Fit astrometry and write the fitted result.",
        dtype=bool,
        default=True
    )
    doPhotometry = pexConfig.Field(
        doc="Fit photometry and write the fitted result.",
        dtype=bool,
        default=True
    )
    coaddName = pexConfig.Field(
        doc="Type of coadd, typically deep or goodSeeing",
        dtype=str,
        default="deep"
    )
    posError = pexConfig.Field(
        doc="Constant term for error on position (in pixel unit)",
        dtype=float,
        default=0.02,
    )
    # TODO: DM-6885 matchCut should be an afw.geom.Angle
    matchCut = pexConfig.Field(
        doc="Matching radius between fitted and reference stars (arcseconds)",
        dtype=float,
        default=3.0,
    )
    minMeasurements = pexConfig.Field(
        doc="Minimum number of associated measured stars for a fitted star to be included in the fit",
        dtype=int,
        default=2,
    )
    minMeasuredStarsPerCcd = pexConfig.Field(
        doc="Minimum number of measuredStars per ccdImage before printing warnings",
        dtype=int,
        default=100,
    )
    minRefStarsPerCcd = pexConfig.Field(
        doc="Minimum number of measuredStars per ccdImage before printing warnings",
        dtype=int,
        default=30,
    )
    astrometrySimpleOrder = pexConfig.Field(
        doc="Polynomial order for fitting the simple astrometry model.",
        dtype=int,
        default=3,
    )
    astrometryChipOrder = pexConfig.Field(
        doc="Order of the per-chip transform for the constrained astrometry model.",
        dtype=int,
        default=1,
    )
    astrometryVisitOrder = pexConfig.Field(
        doc="Order of the per-visit transform for the constrained astrometry model.",
        dtype=int,
        default=5,
    )
    useInputWcs = pexConfig.Field(
        doc="Use the input calexp WCSs to initialize a SimpleAstrometryModel.",
        dtype=bool,
        default=True,
    )
    astrometryModel = pexConfig.ChoiceField(
        doc="Type of model to fit to astrometry",
        dtype=str,
        default="simple",
        allowed={"simple": "One polynomial per ccd",
                 "constrained": "One polynomial per ccd, and one polynomial per visit"}
    )
    photometryModel = pexConfig.ChoiceField(
        doc="Type of model to fit to photometry",
        dtype=str,
        default="simple",
        allowed={"simple": "One constant zeropoint per ccd and visit",
                 "constrained": "Constrained zeropoint per ccd, and one polynomial per visit"}
    )
    photometryVisitOrder = pexConfig.Field(
        doc="Order of the per-visit polynomial transform for the constrained photometry model.",
        dtype=int,
        default=7,
    )
    photometryDoRankUpdate = pexConfig.Field(
        doc="Do the rank update step during minimization. "
        "Skipping this can help deal with models that are too non-linear.",
        dtype=bool,
        default=True,
    )
    astrometryDoRankUpdate = pexConfig.Field(
        doc="Do the rank update step during minimization (should not change the astrometry fit). "
        "Skipping this can help deal with models that are too non-linear.",
        dtype=bool,
        default=True,
    )
    outlierRejectSigma = pexConfig.Field(
        doc="How many sigma to reject outliers at during minimization.",
        dtype=float,
        default=5.0,
    )
    maxPhotometrySteps = pexConfig.Field(
        doc="Maximum number of minimize iterations to take when fitting photometry.",
        dtype=int,
        default=20,
    )
    maxAstrometrySteps = pexConfig.Field(
        doc="Maximum number of minimize iterations to take when fitting photometry.",
        dtype=int,
        default=20,
    )
    astrometryRefObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="Reference object loader for astrometric fit",
    )
    photometryRefObjLoader = pexConfig.ConfigurableField(
        target=LoadIndexedReferenceObjectsTask,
        doc="Reference object loader for photometric fit",
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching",
        default="astrometry"
    )
    writeInitMatrix = pexConfig.Field(
        dtype=bool,
        doc="Write the pre/post-initialization Hessian and gradient to text files, for debugging."
            "The output files will be of the form 'astrometry_preinit-mat.txt', in the current directory."
            "Note that these files are the dense versions of the matrix, and so may be very large.",
        default=False
    )
    writeChi2ContributionFiles = pexConfig.Field(
        dtype=bool,
        doc="Write initial/final fit files containing the contributions to chi2.",
        default=False
    )
    sourceFluxType = pexConfig.Field(
        dtype=str,
        doc="Source flux field to use in source selection and to get fluxes from the catalog.",
        default='Calib'
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["astrometry"]
        sourceSelector.setDefaults()
        # don't want to lose existing flags, just add to them.
        sourceSelector.badFlags.extend(["slot_Shape_flag"])
        # This should be used to set the FluxField value in jointcal::JointcalControl
        sourceSelector.sourceFluxType = self.sourceFluxType


class JointcalTask(pipeBase.CmdLineTask):
    """Jointly astrometrically and photometrically calibrate a group of images."""

    ConfigClass = JointcalConfig
    RunnerClass = JointcalRunner
    _DefaultName = "jointcal"

    def __init__(self, butler=None, profile_jointcal=False, **kwargs):
        """
        Instantiate a JointcalTask.

        Parameters
        ----------
        butler : lsst.daf.persistence.Butler
            The butler is passed to the refObjLoader constructor in case it is
            needed. Ignored if the refObjLoader argument provides a loader directly.
            Used to initialize the astrometry and photometry refObjLoaders.
        profile_jointcal : bool
            set to True to profile different stages of this jointcal run.
        """
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.profile_jointcal = profile_jointcal
        self.makeSubtask("sourceSelector")
        if self.config.doAstrometry:
            self.makeSubtask('astrometryRefObjLoader', butler=butler)
        if self.config.doPhotometry:
            self.makeSubtask('photometryRefObjLoader', butler=butler)

        # To hold various computed metrics for use by tests
        self.job = Job.load_metrics_package(subset='jointcal')

    # We don't need to persist config and metadata at this stage.
    # In this way, we don't need to put a specific entry in the camera mapper policy file
    def _getConfigName(self):
        return None

    def _getMetadataName(self):
        return None

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser"""
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_argument("--profile_jointcal", default=False, action="store_true",
                            help="Profile steps of jointcal separately.")
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=6789 ccd=0..9",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    def _build_ccdImage(self, dataRef, associations, jointcalControl):
        """
        Extract the necessary things from this dataRef to add a new ccdImage.

        Parameters
        ----------
        dataRef : lsst.daf.persistence.ButlerDataRef
            dataRef to extract info from.
        associations : lsst.jointcal.Associations
            object to add the info to, to construct a new CcdImage
        jointcalControl : jointcal.JointcalControl
            control object for associations management

        Returns
        ------
        namedtuple
            wcs : lsst.afw.geom.SkyWcs
                the TAN WCS of this image, read from the calexp
            key : namedtuple
                a key to identify this dataRef by its visit and ccd ids
            filter : str
                this calexp's filter
        """
        if "visit" in dataRef.dataId.keys():
            visit = dataRef.dataId["visit"]
        else:
            visit = dataRef.getButler().queryMetadata("calexp", ("visit"), dataRef.dataId)[0]

        src = dataRef.get("src", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS, immediate=True)

        visitInfo = dataRef.get('calexp_visitInfo')
        detector = dataRef.get('calexp_detector')
        ccdId = detector.getId()
        calib = dataRef.get('calexp_calib')
        tanWcs = dataRef.get('calexp_wcs')
        bbox = dataRef.get('calexp_bbox')
        filt = dataRef.get('calexp_filter')
        filterName = filt.getName()
        fluxMag0 = calib.getFluxMag0()
        photoCalib = afwImage.PhotoCalib(1.0/fluxMag0[0], fluxMag0[1]/fluxMag0[0]**2, bbox)

        goodSrc = self.sourceSelector.run(src)

        if len(goodSrc.sourceCat) == 0:
            self.log.warn("No sources selected in visit %s ccd %s", visit, ccdId)
        else:
            self.log.info("%d sources selected in visit %d ccd %d", len(goodSrc.sourceCat), visit, ccdId)
        associations.createCcdImage(goodSrc.sourceCat,
                                    tanWcs,
                                    visitInfo,
                                    bbox,
                                    filterName,
                                    photoCalib,
                                    detector,
                                    visit,
                                    ccdId,
                                    jointcalControl)

        Result = collections.namedtuple('Result_from_build_CcdImage', ('wcs', 'key', 'filter'))
        Key = collections.namedtuple('Key', ('visit', 'ccd'))
        return Result(tanWcs, Key(visit, ccdId), filterName)

    @pipeBase.timeMethod
    def run(self, dataRefs, profile_jointcal=False):
        """
        Jointly calibrate the astrometry and photometry across a set of images.

        Parameters
        ----------
        dataRefs : list of lsst.daf.persistence.ButlerDataRef
            List of data references to the exposures to be fit.
        profile_jointcal : bool
            Profile the individual steps of jointcal.

        Returns
        -------
        pipe.base.Struct
            struct containing:
            * dataRefs: the provided data references that were fit (with updated WCSs)
            * oldWcsList: the original WCS from each dataRef
            * metrics: dictionary of internally-computed metrics for testing/validation.
        """
        if len(dataRefs) == 0:
            raise ValueError('Need a non-empty list of data references!')

        exitStatus = 0  # exit status for shell

        sourceFluxField = "slot_%sFlux" % (self.config.sourceFluxType,)
        jointcalControl = lsst.jointcal.JointcalControl(sourceFluxField)
        associations = lsst.jointcal.Associations()

        visit_ccd_to_dataRef = {}
        oldWcsList = []
        filters = []
        load_cat_prof_file = 'jointcal_build_ccdImage.prof' if profile_jointcal else ''
        with pipeBase.cmdLineTask.profile(load_cat_prof_file):
            # We need the bounding-box of the focal plane for photometry visit models.
            # NOTE: we only need to read it once, because its the same for all exposures of a camera.
            camera = dataRefs[0].get('camera', immediate=True)
            self.focalPlaneBBox = camera.getFpBBox()
            for ref in dataRefs:
                result = self._build_ccdImage(ref, associations, jointcalControl)
                oldWcsList.append(result.wcs)
                visit_ccd_to_dataRef[result.key] = ref
                filters.append(result.filter)
        filters = collections.Counter(filters)

        associations.computeCommonTangentPoint()

        # Use external reference catalogs handled by LSST stack mechanism
        # Get the bounding box overlapping all associated images
        # ==> This is probably a bad idea to do it this way <== To be improved
        bbox = associations.getRaDecBBox()
        # with Python 3 this can be simplified to afwGeom.SpherePoint(*bbox.getCenter(), afwGeom.degrees)
        bboxCenter = bbox.getCenter()
        center = afwGeom.SpherePoint(bboxCenter[0], bboxCenter[1], afwGeom.degrees)
        bboxMax = bbox.getMax()
        corner = afwGeom.SpherePoint(bboxMax[0], bboxMax[1], afwGeom.degrees)
        radius = center.separation(corner).asRadians()

        # Get astrometry_net_data path
        anDir = lsst.utils.getPackageDir('astrometry_net_data')
        if anDir is None:
            raise RuntimeError("astrometry_net_data is not setup")

        # Determine a default filter associated with the catalog. See DM-9093
        defaultFilter = filters.most_common(1)[0][0]
        self.log.debug("Using %s band for reference flux", defaultFilter)

        # TODO: need a better way to get the tract.
        tract = dataRefs[0].dataId['tract']

        if self.config.doAstrometry:
            astrometry = self._do_load_refcat_and_fit(associations, defaultFilter, center, radius,
                                                      name="astrometry",
                                                      refObjLoader=self.astrometryRefObjLoader,
                                                      fit_function=self._fit_astrometry,
                                                      profile_jointcal=profile_jointcal,
                                                      tract=tract)
            self._write_astrometry_results(associations, astrometry.model, visit_ccd_to_dataRef)
        else:
            astrometry = Astrometry(None, None, None)

        if self.config.doPhotometry:
            photometry = self._do_load_refcat_and_fit(associations, defaultFilter, center, radius,
                                                      name="photometry",
                                                      refObjLoader=self.photometryRefObjLoader,
                                                      fit_function=self._fit_photometry,
                                                      profile_jointcal=profile_jointcal,
                                                      tract=tract,
                                                      filters=filters,
                                                      reject_bad_fluxes=True)
            self._write_photometry_results(associations, photometry.model, visit_ccd_to_dataRef)
        else:
            photometry = Photometry(None, None)

        return pipeBase.Struct(dataRefs=dataRefs,
                               oldWcsList=oldWcsList,
                               job=self.job,
                               exitStatus=exitStatus)

    def _do_load_refcat_and_fit(self, associations, defaultFilter, center, radius,
                                name="", refObjLoader=None, filters=[], fit_function=None,
                                tract=None, profile_jointcal=False, match_cut=3.0,
                                reject_bad_fluxes=False):
        """Load reference catalog, perform the fit, and return the result.

        Parameters
        ----------
        associations : lsst.jointcal.Associations
            The star/reference star associations to fit.
        defaultFilter : str
            filter to load from reference catalog.
        center : lsst.afw.geom.SpherePoint
            ICRS center of field to load from reference catalog.
        radius : lsst.afw.geom.Angle
            On-sky radius to load from reference catalog.
        name : str
            Name of thing being fit: "Astrometry" or "Photometry".
        refObjLoader : lsst.meas.algorithms.LoadReferenceObjectsTask
            Reference object loader to load from for fit.
        filters : list of str, optional
            List of filters to load from the reference catalog.
        fit_function : function
            function to call to perform fit (takes associations object).
        tract : str
            Name of tract currently being fit.
        profile_jointcal : bool, optional
            Separately profile the fitting step.
        match_cut : float, optional
            Radius in arcseconds to find cross-catalog matches to during
            associations.associateCatalogs.
        reject_bad_fluxes : bool, optional
            Reject refCat sources with NaN/inf flux or NaN/0 fluxErr.

        Returns
        -------
        Result of `fit_function()`
        """
        self.log.info("====== Now processing %s...", name)
        # TODO: this should not print "trying to invert a singular transformation:"
        # if it does that, something's not right about the WCS...
        associations.associateCatalogs(match_cut)
        add_measurement(self.job, 'jointcal.associated_%s_fittedStars' % name,
                        associations.fittedStarListSize())

        skyCircle = refObjLoader.loadSkyCircle(center,
                                               afwGeom.Angle(radius, afwGeom.radians),
                                               defaultFilter)

        # Need memory contiguity to get reference filters as a vector.
        if not skyCircle.refCat.isContiguous():
            refCat = skyCircle.refCat.copy(deep=True)
        else:
            refCat = skyCircle.refCat

        # load the reference catalog fluxes.
        # TODO: Simon will file a ticket for making this better (and making it use the color terms)
        refFluxes = {}
        refFluxErrs = {}
        for filt in filters:
            filtKeys = lsst.meas.algorithms.getRefFluxKeys(refCat.schema, filt)
            refFluxes[filt] = refCat.get(filtKeys[0])
            refFluxErrs[filt] = refCat.get(filtKeys[1])

        associations.collectRefStars(refCat, self.config.matchCut*afwGeom.arcseconds,
                                     skyCircle.fluxField, refFluxes, refFluxErrs, reject_bad_fluxes)
        add_measurement(self.job, 'jointcal.collected_%s_refStars' % name,
                        associations.refStarListSize())

        associations.prepareFittedStars(self.config.minMeasurements)

        self._check_star_lists(associations, name)
        add_measurement(self.job, 'jointcal.selected_%s_refStars' % name,
                        associations.nFittedStarsWithAssociatedRefStar())
        add_measurement(self.job, 'jointcal.selected_%s_fittedStars' % name,
                        associations.fittedStarListSize())
        add_measurement(self.job, 'jointcal.selected_%s_ccdImages' % name,
                        associations.nCcdImagesValidForFit())

        load_cat_prof_file = 'jointcal_fit_%s.prof'%name if profile_jointcal else ''
        dataName = "{}_{}".format(tract, defaultFilter)
        with pipeBase.cmdLineTask.profile(load_cat_prof_file):
            result = fit_function(associations, dataName)
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2ContributionFiles:
            baseName = "{}_final_chi2-{}.csv".format(name, dataName)
            result.fit.saveChi2Contributions(baseName)

        return result

    def _check_star_lists(self, associations, name):
        # TODO: these should be len(blah), but we need this properly wrapped first.
        if associations.nCcdImagesValidForFit() == 0:
            raise RuntimeError('No images in the ccdImageList!')
        if associations.fittedStarListSize() == 0:
            raise RuntimeError('No stars in the {} fittedStarList!'.format(name))
        if associations.refStarListSize() == 0:
            raise RuntimeError('No stars in the {} reference star list!'.format(name))

    def _fit_photometry(self, associations, dataName=None):
        """
        Fit the photometric data.

        Parameters
        ----------
        associations : lsst.jointcal.Associations
            The star/reference star associations to fit.
        dataName : str
            Name of the data being processed (e.g. "1234_HSC-Y"), for
            identifying debugging files.

        Returns
        -------
        namedtuple
            fit : lsst.jointcal.PhotometryFit
                The photometric fitter used to perform the fit.
            model : lsst.jointcal.PhotometryModel
                The photometric model that was fit.
        """
        self.log.info("=== Starting photometric fitting...")

        # TODO: should use pex.config.RegistryField here (see DM-9195)
        if self.config.photometryModel == "constrained":
            model = lsst.jointcal.ConstrainedPhotometryModel(associations.getCcdImageList(),
                                                             self.focalPlaneBBox,
                                                             visitOrder=self.config.photometryVisitOrder)
        elif self.config.photometryModel == "simple":
            model = lsst.jointcal.SimpleFluxModel(associations.getCcdImageList())

        fit = lsst.jointcal.PhotometryFit(associations, model)
        chi2 = fit.computeChi2()
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2ContributionFiles:
            baseName = "photometry_initial_chi2-{}.csv".format(dataName)
            fit.saveChi2Contributions(baseName)

        if not np.isfinite(chi2.chi2):
            raise FloatingPointError('Initial chi2 is invalid: %s'%chi2)
        self.log.info("Initialized: %s", str(chi2))
        # The constrained model needs the visit transfo fit first; the chip
        # transfo is initialized from the singleFrame PhotoCalib, so it's close.
        dumpMatrixFile = "photometry_preinit" if self.config.writeInitMatrix else ""
        if self.config.photometryModel == "constrained":
            # TODO: (related to DM-8046): implement Visit/Chip choice
            fit.minimize("ModelVisit", dumpMatrixFile=dumpMatrixFile)
            chi2 = fit.computeChi2()
            self.log.info(str(chi2))
            dumpMatrixFile = ""  # so we don't redo the output on the next step
        fit.minimize("Model", dumpMatrixFile=dumpMatrixFile)
        chi2 = fit.computeChi2()
        self.log.info(str(chi2))
        fit.minimize("Fluxes")
        chi2 = fit.computeChi2()
        self.log.info(str(chi2))
        fit.minimize("Model Fluxes")
        chi2 = fit.computeChi2()
        if not np.isfinite(chi2.chi2):
            raise FloatingPointError('Pre-iteration chi2 is invalid: %s'%chi2)
        self.log.info("Fit prepared with %s", str(chi2))

        model.freezeErrorTransform()
        self.log.debug("Photometry error scales are frozen.")

        chi2 = self._iterate_fit(associations,
                                 fit,
                                 model,
                                 self.config.maxPhotometrySteps,
                                 "photometry",
                                 "Model Fluxes",
                                 doRankUpdate=self.config.photometryDoRankUpdate)

        add_measurement(self.job, 'jointcal.photometry_final_chi2', chi2.chi2)
        add_measurement(self.job, 'jointcal.photometry_final_ndof', chi2.ndof)
        return Photometry(fit, model)

    def _fit_astrometry(self, associations, dataName=None):
        """
        Fit the astrometric data.

        Parameters
        ----------
        associations : lsst.jointcal.Associations
            The star/reference star associations to fit.
        dataName : str
            Name of the data being processed (e.g. "1234_HSC-Y"), for
            identifying debugging files.

        Returns
        -------
        namedtuple
            fit : lsst.jointcal.AstrometryFit
                The astrometric fitter used to perform the fit.
            model : lsst.jointcal.AstrometryModel
                The astrometric model that was fit.
            sky_to_tan_projection : lsst.jointcal.ProjectionHandler
                The model for the sky to tangent plane projection that was used in the fit.
        """

        self.log.info("=== Starting astrometric fitting...")

        associations.deprojectFittedStars()

        # NOTE: need to return sky_to_tan_projection so that it doesn't get garbage collected.
        # TODO: could we package sky_to_tan_projection and model together so we don't have to manage
        # them so carefully?
        sky_to_tan_projection = lsst.jointcal.OneTPPerVisitHandler(associations.getCcdImageList())

        if self.config.astrometryModel == "constrained":
            model = lsst.jointcal.ConstrainedAstrometryModel(associations.getCcdImageList(),
                                                             sky_to_tan_projection,
                                                             chipOrder=self.config.astrometryChipOrder,
                                                             visitOrder=self.config.astrometryVisitOrder)
        elif self.config.astrometryModel == "simple":
            model = lsst.jointcal.SimpleAstrometryModel(associations.getCcdImageList(),
                                                        sky_to_tan_projection,
                                                        self.config.useInputWcs,
                                                        nNotFit=0,
                                                        order=self.config.astrometrySimpleOrder)

        fit = lsst.jointcal.AstrometryFit(associations, model, self.config.posError)
        chi2 = fit.computeChi2()
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2ContributionFiles:
            baseName = "astrometry_initial_chi2-{}.csv".format(dataName)
            fit.saveChi2Contributions(baseName)

        if not np.isfinite(chi2.chi2):
            raise FloatingPointError('Initial chi2 is invalid: %s'%chi2)
        self.log.info("Initialized: %s", str(chi2))
        dumpMatrixFile = "astrometry_preinit" if self.config.writeInitMatrix else ""
        # The constrained model needs the visit transfo fit first; the chip
        # transfo is initialized from the detector's cameraGeom, so it's close.
        if self.config.astrometryModel == "constrained":
            fit.minimize("DistortionsVisit", dumpMatrixFile=dumpMatrixFile)
            chi2 = fit.computeChi2()
            self.log.info(str(chi2))
            dumpMatrixFile = ""  # so we don't redo the output on the next step
        fit.minimize("Distortions", dumpMatrixFile=dumpMatrixFile)
        chi2 = fit.computeChi2()
        self.log.info(str(chi2))
        fit.minimize("Positions")
        chi2 = fit.computeChi2()
        self.log.info(str(chi2))
        fit.minimize("Distortions Positions")
        chi2 = fit.computeChi2()
        self.log.info(str(chi2))
        if not np.isfinite(chi2.chi2):
            raise FloatingPointError('Pre-iteration chi2 is invalid: %s'%chi2)
        self.log.info("Fit prepared with %s", str(chi2))

        chi2 = self._iterate_fit(associations,
                                 fit,
                                 model,
                                 self.config.maxAstrometrySteps,
                                 "astrometry",
                                 "Distortions Positions",
                                 doRankUpdate=self.config.astrometryDoRankUpdate)

        add_measurement(self.job, 'jointcal.astrometry_final_chi2', chi2.chi2)
        add_measurement(self.job, 'jointcal.astrometry_final_ndof', chi2.ndof)

        return Astrometry(fit, model, sky_to_tan_projection)

    def _check_stars(self, associations):
        """Count measured and reference stars per ccd and warn/log them."""
        for ccdImage in associations.getCcdImageList():
            nMeasuredStars, nRefStars = ccdImage.countStars()
            self.log.debug("ccdImage %s has %s measured and %s reference stars",
                           ccdImage.getName(), nMeasuredStars, nRefStars)
            if nMeasuredStars < self.config.minMeasuredStarsPerCcd:
                self.log.warn("ccdImage %s has only %s measuredStars (desired %s)",
                              ccdImage.getName(), nMeasuredStars, self.config.minMeasuredStarsPerCcd)
            if nRefStars < self.config.minRefStarsPerCcd:
                self.log.warn("ccdImage %s has only %s RefStars (desired %s)",
                              ccdImage.getName(), nRefStars, self.config.minRefStarsPerCcd)

    def _iterate_fit(self, associations, fit, model, max_steps, name, whatToFit, doRankUpdate=True):
        """Run fit.minimize up to max_steps times, returning the final chi2."""

        dumpMatrixFile = "%s_postinit" % name if self.config.writeInitMatrix else ""
        for i in range(max_steps):
            # outlier removal at 5 sigma.
            r = fit.minimize(whatToFit,
                             self.config.outlierRejectSigma,
                             doRankUpdate=doRankUpdate,
                             dumpMatrixFile=dumpMatrixFile)
            dumpMatrixFile = ""  # clear it so we don't write the matrix again.
            chi2 = fit.computeChi2()
            self._check_stars(associations)
            if not np.isfinite(chi2.chi2):
                raise FloatingPointError('Fit iteration chi2 is invalid: %s'%chi2)
            self.log.info(str(chi2))
            if r == MinimizeResult.Converged:
                if doRankUpdate:
                    self.log.debug("fit has converged - no more outliers - redo minimization "
                                   "one more time in case we have lost accuracy in rank update.")
                    # Redo minimization one more time in case we have lost accuracy in rank update
                    r = fit.minimize(whatToFit, 5)  # outliers removal at 5 sigma.
                chi2 = fit.computeChi2()
                self.log.info("Fit completed with: %s", str(chi2))
                break
            elif r == MinimizeResult.Chi2Increased:
                self.log.warn("still some ouliers but chi2 increases - retry")
            elif r == MinimizeResult.Failed:
                raise RuntimeError("Chi2 minimization failure, cannot complete fit.")
            else:
                raise RuntimeError("Unxepected return code from minimize().")
        else:
            self.log.error("%s failed to converge after %d steps"%(name, max_steps))

        return chi2

    def _write_astrometry_results(self, associations, model, visit_ccd_to_dataRef):
        """
        Write the fitted astrometric results to a new 'jointcal_wcs' dataRef.

        Parameters
        ----------
        associations : lsst.jointcal.Associations
            The star/reference star associations to fit.
        model : lsst.jointcal.AstrometryModel
            The astrometric model that was fit.
        visit_ccd_to_dataRef : dict of Key: lsst.daf.persistence.ButlerDataRef
            dict of ccdImage identifiers to dataRefs that were fit
        """

        ccdImageList = associations.getCcdImageList()
        for ccdImage in ccdImageList:
            # TODO: there must be a better way to identify this ccdImage than a visit,ccd pair?
            ccd = ccdImage.ccdId
            visit = ccdImage.visit
            dataRef = visit_ccd_to_dataRef[(visit, ccd)]
            self.log.info("Updating WCS for visit: %d, ccd: %d", visit, ccd)
            skyWcs = model.makeSkyWcs(ccdImage)
            try:
                dataRef.put(skyWcs, 'jointcal_wcs')
            except pexExceptions.Exception as e:
                self.log.fatal('Failed to write updated Wcs: %s', str(e))
                raise e

    def _write_photometry_results(self, associations, model, visit_ccd_to_dataRef):
        """
        Write the fitted photometric results to a new 'jointcal_photoCalib' dataRef.

        Parameters
        ----------
        associations : lsst.jointcal.Associations
            The star/reference star associations to fit.
        model : lsst.jointcal.PhotometryModel
            The photoometric model that was fit.
        visit_ccd_to_dataRef : dict of Key: lsst.daf.persistence.ButlerDataRef
            dict of ccdImage identifiers to dataRefs that were fit
        """

        ccdImageList = associations.getCcdImageList()
        for ccdImage in ccdImageList:
            # TODO: there must be a better way to identify this ccdImage than a visit,ccd pair?
            ccd = ccdImage.ccdId
            visit = ccdImage.visit
            dataRef = visit_ccd_to_dataRef[(visit, ccd)]
            self.log.info("Updating PhotoCalib for visit: %d, ccd: %d", visit, ccd)
            photoCalib = model.toPhotoCalib(ccdImage)
            try:
                dataRef.put(photoCalib, 'jointcal_photoCalib')
            except pexExceptions.Exception as e:
                self.log.fatal('Failed to write updated PhotoCalib: %s', str(e))
                raise e
