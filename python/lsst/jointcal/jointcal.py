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

import dataclasses
import collections
import os
import logging

import astropy.time
import numpy as np
import astropy.units as u

import lsst.geom
import lsst.utils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.afw.image import fluxErrFromABMagErr
import lsst.afw.cameraGeom
import lsst.afw.table
from lsst.pipe.tasks.colorterms import ColortermLibrary
from lsst.verify import Job, Measurement

from lsst.meas.algorithms import (ReferenceObjectLoader, ReferenceSourceSelectorTask,
                                  LoadReferenceObjectsConfig)
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

import lsst.jointcal
from lsst.jointcal import MinimizeResult

__all__ = ["JointcalConfig", "JointcalTask"]

Photometry = collections.namedtuple('Photometry', ('fit', 'model'))
Astrometry = collections.namedtuple('Astrometry', ('fit', 'model', 'sky_to_tan_projection'))


# TODO: move this to MeasurementSet in lsst.verify per DM-12655.
def add_measurement(job, name, value):
    meas = Measurement(job.metrics[name], value)
    job.measurements.insert(meas)


class JointcalTaskConnections(pipeBase.PipelineTaskConnections,
                              dimensions=("skymap", "tract", "instrument", "physical_filter")):
    """Middleware input/output connections for jointcal data."""
    inputCamera = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The camera instrument that took these observations.",
        name="camera",
        storageClass="Camera",
        dimensions=("instrument",),
        isCalibration=True,
    )
    inputSourceTableVisit = pipeBase.connectionTypes.Input(
        doc="Source table in parquet format, per visit",
        name="sourceTable_visit",
        storageClass="DataFrame",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )
    inputVisitSummary = pipeBase.connectionTypes.Input(
        doc=("Per-visit consolidated exposure metadata built from calexps. "
             "These catalogs use detector id for the id and must be sorted for "
             "fast lookups of a detector."),
        name="visitSummary",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
        deferLoad=True,
        multiple=True,
    )
    astrometryRefCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The astrometry reference catalog to match to loaded input catalog sources.",
        name="gaia_dr2_20200414",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )
    photometryRefCat = pipeBase.connectionTypes.PrerequisiteInput(
        doc="The photometry reference catalog to match to loaded input catalog sources.",
        name="ps1_pv3_3pi_20170110",
        storageClass="SimpleCatalog",
        dimensions=("skypix",),
        deferLoad=True,
        multiple=True,
    )

    outputWcs = pipeBase.connectionTypes.Output(
        doc=("Per-tract, per-visit world coordinate systems derived from the fitted model."
             " These catalogs only contain entries for detectors with an output, and use"
             " the detector id for the catalog id, sorted on id for fast lookups of a detector."),
        name="jointcalSkyWcsCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit", "skymap", "tract"),
        multiple=True
    )
    outputPhotoCalib = pipeBase.connectionTypes.Output(
        doc=("Per-tract, per-visit photometric calibrations derived from the fitted model."
             " These catalogs only contain entries for detectors with an output, and use"
             " the detector id for the catalog id, sorted on id for fast lookups of a detector."),
        name="jointcalPhotoCalibCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit", "skymap", "tract"),
        multiple=True
    )

    # measurements of metrics
    # The vars() trick used here allows us to set class attributes
    # programatically. Taken from:
    # https://stackoverflow.com/questions/2519807/setting-a-class-attribute-with-a-given-name-in-python-while-defining-the-class
    for name in ("astrometry", "photometry"):
        vars()[f"{name}_matched_fittedStars"] = pipeBase.connectionTypes.Output(
            doc=f"The number of cross-matched fittedStars for {name}",
            name=f"metricvalue_jointcal_{name}_matched_fittedStars",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_collected_refStars"] = pipeBase.connectionTypes.Output(
            doc=f"The number of {name} reference stars drawn from the reference catalog, before matching.",
            name=f"metricvalue_jointcal_{name}_collected_refStars",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_prepared_refStars"] = pipeBase.connectionTypes.Output(
            doc=f"The number of {name} reference stars matched to fittedStars.",
            name=f"metricvalue_jointcal_{name}_prepared_refStars",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_prepared_fittedStars"] = pipeBase.connectionTypes.Output(
            doc=f"The number of cross-matched fittedStars after cleanup, for {name}.",
            name=f"metricvalue_jointcal_{name}_prepared_fittedStars",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_prepared_ccdImages"] = pipeBase.connectionTypes.Output(
            doc=f"The number of ccdImages that will be fit for {name}, after cleanup.",
            name=f"metricvalue_jointcal_{name}_prepared_ccdImages",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_final_chi2"] = pipeBase.connectionTypes.Output(
            doc=f"The final chi2 of the {name} fit.",
            name=f"metricvalue_jointcal_{name}_final_chi2",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )
        vars()[f"{name}_final_ndof"] = pipeBase.connectionTypes.Output(
            doc=f"The number of degrees of freedom of the fitted {name}.",
            name=f"metricvalue_jointcal_{name}_final_ndof",
            storageClass="MetricValue",
            dimensions=("skymap", "tract", "instrument", "physical_filter"),
        )

    def __init__(self, *, config=None):
        super().__init__(config=config)
        # When we are only doing one of astrometry or photometry, we don't
        # need the reference catalog or produce the outputs for the other.
        # This informs the middleware of that when the QuantumGraph is
        # generated, so we don't block on getting something we won't need or
        # create an expectation that downstream tasks will be able to consume
        # something we won't produce.
        if not config.doAstrometry:
            self.prerequisiteInputs.remove("astrometryRefCat")
            self.outputs.remove("outputWcs")
            for key in list(self.outputs):
                if "metricvalue_jointcal_astrometry" in key:
                    self.outputs.remove(key)
        if not config.doPhotometry:
            self.prerequisiteInputs.remove("photometryRefCat")
            self.outputs.remove("outputPhotoCalib")
            for key in list(self.outputs):
                if "metricvalue_jointcal_photometry" in key:
                    self.outputs.remove(key)

    def getSpatialBoundsConnections(self):
        return ("inputVisitSummary",)


class JointcalConfig(pipeBase.PipelineTaskConfig,
                     pipelineConnections=JointcalTaskConnections):
    """Configuration for JointcalTask"""

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
    sourceFluxType = pexConfig.Field(
        dtype=str,
        doc="Source flux field to use in source selection and to get fluxes from the catalog.",
        default='apFlux_12_0'
    )
    positionErrorPedestal = pexConfig.Field(
        doc="Systematic term to apply to the measured position error (pixels)",
        dtype=float,
        default=0.02,
    )
    photometryErrorPedestal = pexConfig.Field(
        doc="Systematic term to apply to the measured error on flux or magnitude as a "
        "fraction of source flux or magnitude delta (e.g. 0.05 is 5% of flux or +50 millimag).",
        dtype=float,
        default=0.0,
    )
    # TODO: DM-6885 matchCut should be an geom.Angle
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
    allowLineSearch = pexConfig.Field(
        doc="Allow a line search during minimization, if it is reasonable for the model"
        " (models with a significant non-linear component, e.g. constrainedPhotometry).",
        dtype=bool,
        default=False
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
        default="constrained",
        allowed={"simple": "One polynomial per ccd",
                 "constrained": "One polynomial per ccd, and one polynomial per visit"}
    )
    photometryModel = pexConfig.ChoiceField(
        doc="Type of model to fit to photometry",
        dtype=str,
        default="constrainedMagnitude",
        allowed={"simpleFlux": "One constant zeropoint per ccd and visit, fitting in flux space.",
                 "constrainedFlux": "Constrained zeropoint per ccd, and one polynomial per visit,"
                 " fitting in flux space.",
                 "simpleMagnitude": "One constant zeropoint per ccd and visit,"
                 " fitting in magnitude space.",
                 "constrainedMagnitude": "Constrained zeropoint per ccd, and one polynomial per visit,"
                 " fitting in magnitude space.",
                 }
    )
    applyColorTerms = pexConfig.Field(
        doc="Apply photometric color terms to reference stars?"
            "Requires that colorterms be set to a ColortermLibrary",
        dtype=bool,
        default=False
    )
    colorterms = pexConfig.ConfigField(
        doc="Library of photometric reference catalog name to color term dict.",
        dtype=ColortermLibrary,
    )
    photometryVisitOrder = pexConfig.Field(
        doc="Order of the per-visit polynomial transform for the constrained photometry model.",
        dtype=int,
        default=7,
    )
    photometryDoRankUpdate = pexConfig.Field(
        doc=("Do the rank update step during minimization. "
             "Skipping this can help deal with models that are too non-linear."),
        dtype=bool,
        default=True,
    )
    astrometryDoRankUpdate = pexConfig.Field(
        doc=("Do the rank update step during minimization (should not change the astrometry fit). "
             "Skipping this can help deal with models that are too non-linear."),
        dtype=bool,
        default=True,
    )
    outlierRejectSigma = pexConfig.Field(
        doc="How many sigma to reject outliers at during minimization.",
        dtype=float,
        default=5.0,
    )
    astrometryOutlierRelativeTolerance = pexConfig.Field(
        doc=("Convergence tolerance for outlier rejection threshold when fitting astrometry. Iterations will "
             "stop when the fractional change in the chi2 cut level is below this value. If tolerance is set "
             "to zero, iterations will continue until there are no more outliers. We suggest a value of 0.002"
             "as a balance between a shorter minimization runtime and achieving a final fitted model that is"
             "close to the solution found when removing all outliers."),
        dtype=float,
        default=0,
    )
    maxPhotometrySteps = pexConfig.Field(
        doc="Maximum number of minimize iterations to take when fitting photometry.",
        dtype=int,
        default=20,
    )
    maxAstrometrySteps = pexConfig.Field(
        doc="Maximum number of minimize iterations to take when fitting astrometry.",
        dtype=int,
        default=20,
    )
    astrometryRefObjLoader = pexConfig.ConfigField(
        dtype=LoadReferenceObjectsConfig,
        doc="Reference object loader for astrometric fit",
    )
    photometryRefObjLoader = pexConfig.ConfigField(
        dtype=LoadReferenceObjectsConfig,
        doc="Reference object loader for photometric fit",
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching",
        default="science"
    )
    astrometryReferenceSelector = pexConfig.ConfigurableField(
        target=ReferenceSourceSelectorTask,
        doc="How to down-select the loaded astrometry reference catalog.",
    )
    photometryReferenceSelector = pexConfig.ConfigurableField(
        target=ReferenceSourceSelectorTask,
        doc="How to down-select the loaded photometry reference catalog.",
    )
    astrometryReferenceErr = pexConfig.Field(
        doc=("Uncertainty on reference catalog coordinates [mas] to use in place of the `coord_*Err` fields. "
             "If None, then raise an exception if the reference catalog is missing coordinate errors. "
             "If specified, overrides any existing `coord_*Err` values."),
        dtype=float,
        default=None,
        optional=True
    )

    # configs for outputting debug information
    writeInitMatrix = pexConfig.Field(
        dtype=bool,
        doc=("Write the pre/post-initialization Hessian and gradient to text files, for debugging. "
             "Output files will be written to `config.debugOutputPath` and will "
             "be of the form 'astrometry_[pre|post]init-TRACT-FILTER-mat.txt'. "
             "Note that these files are the dense versions of the matrix, and so may be very large."),
        default=False
    )
    writeChi2FilesInitialFinal = pexConfig.Field(
        dtype=bool,
        doc=("Write .csv files containing the contributions to chi2 for the initialization and final fit. "
             "Output files will be written to `config.debugOutputPath` and will "
             "be of the form `astrometry_[initial|final]_chi2-TRACT-FILTER."),
        default=False
    )
    writeChi2FilesOuterLoop = pexConfig.Field(
        dtype=bool,
        doc=("Write .csv files containing the contributions to chi2 for the outer fit loop. "
             "Output files will be written to `config.debugOutputPath` and will "
             "be of the form `astrometry_init-NN_chi2-TRACT-FILTER`."),
        default=False
    )
    writeInitialModel = pexConfig.Field(
        dtype=bool,
        doc=("Write the pre-initialization model to text files, for debugging. "
             "Output files will be written to `config.debugOutputPath` and will be "
             "of the form `initial_astrometry_model-TRACT_FILTER.txt`."),
        default=False
    )
    debugOutputPath = pexConfig.Field(
        dtype=str,
        default=".",
        doc=("Path to write debug output files to. Used by "
             "`writeInitialModel`, `writeChi2Files*`, `writeInitMatrix`.")
    )
    detailedProfile = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Output separate profiling information for different parts of jointcal, e.g. data read, fitting"
    )

    def validate(self):
        super().validate()
        if self.doPhotometry and self.applyColorTerms and len(self.colorterms.data) == 0:
            msg = "applyColorTerms=True requires the `colorterms` field be set to a ColortermLibrary."
            raise pexConfig.FieldValidationError(JointcalConfig.colorterms, self, msg)
        if self.doAstrometry and not self.doPhotometry and self.applyColorTerms:
            msg = ("Only doing astrometry, but Colorterms are not applied for astrometry;"
                   "applyColorTerms=True will be ignored.")
            logging.getLogger("lsst.jointcal").warning(msg)

    def setDefaults(self):
        # Use only primary stars.
        self.sourceSelector["science"].doRequirePrimary = True
        # Use only stars because aperture fluxes of galaxies are biased and depend on seeing.
        self.sourceSelector["science"].doUnresolved = True
        self.sourceSelector["science"].unresolved.name = "sizeExtendedness"
        # with dependable signal to noise ratio.
        self.sourceSelector["science"].doSignalToNoise = True
        # Min SNR must be > 0 because jointcal cannot handle negative fluxes,
        # and S/N > 10 to use sources that are not too faint, and thus better measured.
        self.sourceSelector["science"].signalToNoise.minimum = 10.
        # Base SNR selection on `sourceFluxType` because that is the flux that jointcal fits.
        self.sourceSelector["science"].signalToNoise.fluxField = f"{self.sourceFluxType}_instFlux"
        self.sourceSelector["science"].signalToNoise.errField = f"{self.sourceFluxType}_instFluxErr"
        # Do not trust blended sources" aperture fluxes which also depend on seeing.
        self.sourceSelector["science"].doIsolated = True
        self.sourceSelector["science"].isolated.parentName = "parentSourceId"
        self.sourceSelector["science"].isolated.nChildName = "deblend_nChild"
        # Do not trust either flux or centroid measurements with flags,
        # chosen from the usual QA flags for stars)
        self.sourceSelector["science"].doFlags = True
        badFlags = ["pixelFlags_edge",
                    "pixelFlags_saturated",
                    "pixelFlags_interpolatedCenter",
                    "pixelFlags_interpolated",
                    "pixelFlags_crCenter",
                    "pixelFlags_bad",
                    "hsmPsfMoments_flag",
                    f"{self.sourceFluxType}_flag",
                    ]
        self.sourceSelector["science"].flags.bad = badFlags
        self.sourceSelector["science"].doRequireFiniteRaDec = True
        self.sourceSelector["science"].requireFiniteRaDec.raColName = "ra"
        self.sourceSelector["science"].requireFiniteRaDec.decColName = "dec"

        # Use Gaia-DR2 with proper motions for astrometry; phot_g_mean is the
        # primary Gaia band, but is not like any normal photometric band.
        self.astrometryRefObjLoader.requireProperMotion = True
        self.astrometryRefObjLoader.anyFilterMapsToThis = "phot_g_mean"


def writeModel(model, filename, log):
    """Write model to outfile."""
    with open(filename, "w") as file:
        file.write(repr(model))
    log.info("Wrote %s to file: %s", model, filename)


@dataclasses.dataclass
class JointcalInputData:
    """The input data jointcal needs for each detector/visit."""
    visit: int
    """The visit identifier of this exposure."""
    catalog: lsst.afw.table.SourceCatalog
    """The catalog derived from this exposure."""
    visitInfo: lsst.afw.image.VisitInfo
    """The VisitInfo of this exposure."""
    detector: lsst.afw.cameraGeom.Detector
    """The detector of this exposure."""
    photoCalib: lsst.afw.image.PhotoCalib
    """The photometric calibration of this exposure."""
    wcs: lsst.afw.geom.skyWcs
    """The WCS of this exposure."""
    bbox: lsst.geom.Box2I
    """The bounding box of this exposure."""
    filter: lsst.afw.image.FilterLabel
    """The filter of this exposure."""


class JointcalTask(pipeBase.PipelineTask):
    """Astrometricly and photometricly calibrate across multiple visits of the
    same field.
    """

    ConfigClass = JointcalConfig
    _DefaultName = "jointcal"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask("sourceSelector")
        if self.config.doAstrometry:
            self.makeSubtask("astrometryReferenceSelector")
        else:
            self.astrometryRefObjLoader = None
        if self.config.doPhotometry:
            self.makeSubtask("photometryReferenceSelector")
        else:
            self.photometryRefObjLoader = None

        # To hold various computed metrics for use by tests
        self.job = Job.load_metrics_package(subset='jointcal')

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        # We override runQuantum to set up the refObjLoaders and write the
        # outputs to the correct refs.
        inputs = butlerQC.get(inputRefs)
        # We want the tract number for writing debug files
        tract = butlerQC.quantum.dataId['tract']
        if self.config.doAstrometry:
            self.astrometryRefObjLoader = ReferenceObjectLoader(
                dataIds=[ref.datasetRef.dataId
                         for ref in inputRefs.astrometryRefCat],
                refCats=inputs.pop('astrometryRefCat'),
                config=self.config.astrometryRefObjLoader,
                name=self.config.connections.astrometryRefCat,
                log=self.log)
        if self.config.doPhotometry:
            self.photometryRefObjLoader = ReferenceObjectLoader(
                dataIds=[ref.datasetRef.dataId
                         for ref in inputRefs.photometryRefCat],
                refCats=inputs.pop('photometryRefCat'),
                config=self.config.photometryRefObjLoader,
                name=self.config.connections.photometryRefCat,
                log=self.log)
        outputs = self.run(**inputs, tract=tract)
        self._put_metrics(butlerQC, outputs.job, outputRefs)
        if self.config.doAstrometry:
            self._put_output(butlerQC, outputs.outputWcs, outputRefs.outputWcs,
                             inputs['inputCamera'], "setWcs")
        if self.config.doPhotometry:
            self._put_output(butlerQC, outputs.outputPhotoCalib, outputRefs.outputPhotoCalib,
                             inputs['inputCamera'], "setPhotoCalib")

    def _put_metrics(self, butlerQC, job, outputRefs):
        """Persist all measured metrics stored in a job.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`; This is the input to `runQuantum`.
        job : `lsst.verify.job.Job`
            Measurements of metrics to persist.
        outputRefs : `list` [`lsst.pipe.base.OutputQuantizedConnection`]
            The DatasetRefs to persist the data to.
        """
        for key in job.measurements.keys():
            butlerQC.put(job.measurements[key], getattr(outputRefs, key.fqn.replace('jointcal.', '')))

    def _put_output(self, butlerQC, outputs, outputRefs, camera, setter):
        """Persist the output datasets to their appropriate datarefs.

        Parameters
        ----------
        butlerQC : `lsst.pipe.base.QuantumContext`
            A butler which is specialized to operate in the context of a
            `lsst.daf.butler.Quantum`; This is the input to `runQuantum`.
        outputs : `dict` [`tuple`, `lsst.afw.geom.SkyWcs`] or
                  `dict` [`tuple, `lsst.afw.image.PhotoCalib`]
            The fitted objects to persist.
        outputRefs : `list` [`lsst.pipe.base.OutputQuantizedConnection`]
            The DatasetRefs to persist the data to.
        camera : `lsst.afw.cameraGeom.Camera`
            The camera for this instrument, to get detector ids from.
        setter : `str`
            The method to call on the ExposureCatalog to set each output.
        """
        schema = lsst.afw.table.ExposureTable.makeMinimalSchema()
        schema.addField('visit', type='L', doc='Visit number')

        def new_catalog(visit, size):
            """Return an catalog ready to be filled with appropriate output."""
            catalog = lsst.afw.table.ExposureCatalog(schema)
            catalog.resize(size)
            catalog['visit'] = visit
            metadata = lsst.daf.base.PropertyList()
            metadata.add("COMMENT", "Catalog id is detector id, sorted.")
            metadata.add("COMMENT", "Only detectors with data have entries.")
            return catalog

        # count how many detectors have output for each visit
        detectors_per_visit = collections.defaultdict(int)
        for key in outputs:
            # key is (visit, detector_id), and we only need visit here
            detectors_per_visit[key[0]] += 1

        for ref in outputRefs:
            visit = ref.dataId['visit']
            catalog = new_catalog(visit, detectors_per_visit[visit])
            # Iterate over every detector and skip the ones we don't have output for.
            i = 0
            for detector in camera:
                detectorId = detector.getId()
                key = (ref.dataId['visit'], detectorId)
                if key not in outputs:
                    # skip detectors we don't have output for
                    self.log.debug("No %s output for detector %s in visit %s",
                                   setter[3:], detectorId, visit)
                    continue

                catalog[i].setId(detectorId)
                getattr(catalog[i], setter)(outputs[key])
                i += 1

            catalog.sort()  # ensure that the detectors are in sorted order, for fast lookups
            butlerQC.put(catalog, ref)
            self.log.info("Wrote %s detectors to %s", i, ref)

    def run(self, inputSourceTableVisit, inputVisitSummary, inputCamera, tract=None):
        # Docstring inherited.

        # We take values out of the Parquet table, and put them in "flux_",
        # and the config.sourceFluxType field is used during that extraction,
        # so just use "flux" here.
        sourceFluxField = "flux"
        jointcalControl = lsst.jointcal.JointcalControl(sourceFluxField)
        associations = lsst.jointcal.Associations()
        self.focalPlaneBBox = inputCamera.getFpBBox()
        oldWcsList, bands = self._load_data(inputSourceTableVisit,
                                            inputVisitSummary,
                                            associations,
                                            jointcalControl,
                                            inputCamera)

        boundingCircle, center, radius, defaultFilter, epoch = self._prep_sky(associations, bands)

        if self.config.doAstrometry:
            astrometry = self._do_load_refcat_and_fit(associations, defaultFilter, center, radius,
                                                      name="astrometry",
                                                      refObjLoader=self.astrometryRefObjLoader,
                                                      referenceSelector=self.astrometryReferenceSelector,
                                                      fit_function=self._fit_astrometry,
                                                      tract=tract,
                                                      epoch=epoch)
            astrometry_output = self._make_output(associations.getCcdImageList(),
                                                  astrometry.model,
                                                  "makeSkyWcs")
        else:
            astrometry_output = None

        if self.config.doPhotometry:
            photometry = self._do_load_refcat_and_fit(associations, defaultFilter, center, radius,
                                                      name="photometry",
                                                      refObjLoader=self.photometryRefObjLoader,
                                                      referenceSelector=self.photometryReferenceSelector,
                                                      fit_function=self._fit_photometry,
                                                      tract=tract,
                                                      epoch=epoch,
                                                      reject_bad_fluxes=True)
            photometry_output = self._make_output(associations.getCcdImageList(),
                                                  photometry.model,
                                                  "toPhotoCalib")
        else:
            photometry_output = None

        return pipeBase.Struct(outputWcs=astrometry_output,
                               outputPhotoCalib=photometry_output,
                               job=self.job,
                               astrometryRefObjLoader=self.astrometryRefObjLoader,
                               photometryRefObjLoader=self.photometryRefObjLoader)

    def _load_data(self, inputSourceTableVisit, inputVisitSummary, associations,
                   jointcalControl, camera):
        """Read the data that jointcal needs to run.

        Modifies ``associations`` in-place with the loaded data.

        Parameters
        ----------
        inputSourceTableVisit : `list` [`lsst.daf.butler.DeferredDatasetHandle`]
            References to visit-level DataFrames to load the catalog data from.
        inputVisitSummary : `list` [`lsst.daf.butler.DeferredDatasetHandle`]
            Visit-level exposure summary catalog with metadata.
        associations : `lsst.jointcal.Associations`
            Object to add the loaded data to by constructing new CcdImages.
        jointcalControl : `jointcal.JointcalControl`
            Control object for C++ associations management.
        camera : `lsst.afw.cameraGeom.Camera`
            Camera object for detector geometry.

        Returns
        -------
        oldWcsList: `list` [`lsst.afw.geom.SkyWcs`]
            The original WCS of the input data, to aid in writing tests.
        bands : `list` [`str`]
            The filter bands of each input dataset.
        """
        oldWcsList = []
        filters = []
        load_cat_profile_file = 'jointcal_load_data.prof' if self.config.detailedProfile else ''
        with lsst.utils.timer.profile(load_cat_profile_file):
            table = make_schema_table()  # every detector catalog has the same layout
            # No guarantee that the input is in the same order of visits, so we have to map one of them.
            catalogMap = {ref.dataId['visit']: i for i, ref in enumerate(inputSourceTableVisit)}
            detectorDict = {detector.getId(): detector for detector in camera}

            columns = None

            for visitSummaryRef in inputVisitSummary:
                visitSummary = visitSummaryRef.get()

                dataRef = inputSourceTableVisit[catalogMap[visitSummaryRef.dataId['visit']]]
                if columns is None:
                    inColumns = dataRef.get(component='columns')
                    columns, ixxColumns = get_sourceTable_visit_columns(inColumns,
                                                                        self.config,
                                                                        self.sourceSelector)
                visitCatalog = dataRef.get(parameters={'columns': columns})

                selected = self.sourceSelector.run(visitCatalog)
                if len(selected.sourceCat) == 0:
                    self.log.warning("No sources selected in visit %s.  Skipping...",
                                     visitSummary["visit"][0])
                    continue

                # Build a CcdImage for each detector in this visit.
                detectors = {id: index for index, id in enumerate(visitSummary['id'])}
                for id, index in detectors.items():
                    catalog = extract_detector_catalog_from_visit_catalog(table,
                                                                          selected.sourceCat,
                                                                          id,
                                                                          ixxColumns,
                                                                          self.config.sourceFluxType,
                                                                          self.log)
                    if catalog is None:
                        continue
                    data = self._make_one_input_data(visitSummary[index], catalog, detectorDict)
                    result = self._build_ccdImage(data, associations, jointcalControl)
                    if result is not None:
                        oldWcsList.append(result.wcs)
                        # A visit has only one band, so we can just use the first.
                        filters.append(data.filter)
        if len(filters) == 0:
            raise RuntimeError("No data to process: did source selector remove all sources?")
        filters = collections.Counter(filters)

        return oldWcsList, filters

    def _make_one_input_data(self, visitRecord, catalog, detectorDict):
        """Return a data structure for this detector+visit."""
        return JointcalInputData(visit=visitRecord['visit'],
                                 catalog=catalog,
                                 visitInfo=visitRecord.getVisitInfo(),
                                 detector=detectorDict[visitRecord.getId()],
                                 photoCalib=visitRecord.getPhotoCalib(),
                                 wcs=visitRecord.getWcs(),
                                 bbox=visitRecord.getBBox(),
                                 # ExposureRecord doesn't have a FilterLabel yet,
                                 # so we have to make one.
                                 filter=lsst.afw.image.FilterLabel(band=visitRecord['band'],
                                                                   physical=visitRecord['physical_filter']))

    def _build_ccdImage(self, data, associations, jointcalControl):
        """
        Extract the necessary things from this catalog+metadata to add a new
        ccdImage.

        Parameters
        ----------
        data : `JointcalInputData`
            The loaded input data.
        associations : `lsst.jointcal.Associations`
            Object to add the info to, to construct a new CcdImage
        jointcalControl : `jointcal.JointcalControl`
            Control object for associations management

        Returns
        ------
        namedtuple or `None`
            ``wcs``
                The TAN WCS of this image, read from the calexp
                (`lsst.afw.geom.SkyWcs`).
            ``key``
                A key to identify this dataRef by its visit and ccd ids
                (`namedtuple`).
            `None`
            if there are no sources in the loaded catalog.
        """
        if len(data.catalog) == 0:
            self.log.warning("No sources selected in visit %s ccd %s", data.visit, data.detector.getId())
            return None

        associations.createCcdImage(data.catalog,
                                    data.wcs,
                                    data.visitInfo,
                                    data.bbox,
                                    data.filter.physicalLabel,
                                    data.photoCalib,
                                    data.detector,
                                    data.visit,
                                    data.detector.getId(),
                                    jointcalControl)

        Result = collections.namedtuple('Result_from_build_CcdImage', ('wcs', 'key'))
        Key = collections.namedtuple('Key', ('visit', 'ccd'))
        return Result(data.wcs, Key(data.visit, data.detector.getId()))

    def _getDebugPath(self, filename):
        """Constructs a path to filename using the configured debug path.
        """
        return os.path.join(self.config.debugOutputPath, filename)

    def _prep_sky(self, associations, filters):
        """Prepare on-sky and other data that must be computed after data has
        been read.
        """
        associations.computeCommonTangentPoint()

        boundingCircle = associations.computeBoundingCircle()
        center = lsst.geom.SpherePoint(boundingCircle.getCenter())
        radius = lsst.geom.Angle(boundingCircle.getOpeningAngle().asRadians(), lsst.geom.radians)

        self.log.info(f"Data has center={center} with radius={radius.asDegrees()} degrees.")

        # Determine a default filter band associated with the catalog. See DM-9093
        defaultFilter = filters.most_common(1)[0][0]
        self.log.debug("Using '%s' filter for reference flux", defaultFilter.physicalLabel)

        # compute and set the reference epoch of the observations, for proper motion corrections
        epoch = self._compute_proper_motion_epoch(associations.getCcdImageList())
        associations.setEpoch(epoch.jyear)

        return boundingCircle, center, radius, defaultFilter, epoch

    def _get_refcat_coordinate_error_override(self, refCat, name):
        """Check whether we should override the refcat coordinate errors, and
        return the overridden error if necessary.

        Parameters
        ----------
        refCat : `lsst.afw.table.SimpleCatalog`
            The reference catalog to check for a ``coord_raErr`` field.
        name : `str`
            Whether we are doing "astrometry" or "photometry".

        Returns
        -------
        refCoordErr : `float`
            The refcat coordinate error to use, or NaN if we are not overriding
            those fields.

        Raises
        ------
        lsst.pex.config.FieldValidationError
            Raised if the refcat does not contain coordinate errors and
            ``config.astrometryReferenceErr`` is not set.
        """
        # This value doesn't matter for photometry, so just set something to
        # keep old refcats from causing problems.
        if name.lower() == "photometry":
            if 'coord_raErr' not in refCat.schema:
                return 100
            else:
                return float('nan')

        if self.config.astrometryReferenceErr is None and 'coord_raErr' not in refCat.schema:
            msg = ("Reference catalog does not contain coordinate errors, "
                   "and config.astrometryReferenceErr not supplied.")
            raise pexConfig.FieldValidationError(JointcalConfig.astrometryReferenceErr,
                                                 self.config,
                                                 msg)

        if self.config.astrometryReferenceErr is not None and 'coord_raErr' in refCat.schema:
            self.log.warning("Overriding reference catalog coordinate errors with %f/coordinate [mas]",
                             self.config.astrometryReferenceErr)

        if self.config.astrometryReferenceErr is None:
            return float('nan')
        else:
            return self.config.astrometryReferenceErr

    def _compute_proper_motion_epoch(self, ccdImageList):
        """Return the proper motion correction epoch of the provided images.

        Parameters
        ----------
        ccdImageList : `list` [`lsst.jointcal.CcdImage`]
            The images to compute the appropriate epoch for.

        Returns
        -------
        epoch : `astropy.time.Time`
            The date to use for proper motion corrections.
        """
        return astropy.time.Time(np.mean([ccdImage.getEpoch() for ccdImage in ccdImageList]),
                                 format="jyear",
                                 scale="tai")

    def _do_load_refcat_and_fit(self, associations, defaultFilter, center, radius,
                                tract="", match_cut=3.0,
                                reject_bad_fluxes=False, *,
                                name="", refObjLoader=None, referenceSelector=None,
                                fit_function=None, epoch=None):
        """Load reference catalog, perform the fit, and return the result.

        Parameters
        ----------
        associations : `lsst.jointcal.Associations`
            The star/reference star associations to fit.
        defaultFilter : `lsst.afw.image.FilterLabel`
            filter to load from reference catalog.
        center : `lsst.geom.SpherePoint`
            ICRS center of field to load from reference catalog.
        radius : `lsst.geom.Angle`
            On-sky radius to load from reference catalog.
        name : `str`
            Name of thing being fit: "astrometry" or "photometry".
        refObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader`
            Reference object loader to use to load a reference catalog.
        referenceSelector : `lsst.meas.algorithms.ReferenceSourceSelectorTask`
            Selector to use to pick objects from the loaded reference catalog.
        fit_function : callable
            Function to call to perform fit (takes Associations object).
        tract : `str`, optional
            Name of tract currently being fit.
        match_cut : `float`, optional
            Radius in arcseconds to find cross-catalog matches to during
            associations.associateCatalogs.
        reject_bad_fluxes : `bool`, optional
            Reject refCat sources with NaN/inf flux or NaN/0 fluxErr.
        epoch : `astropy.time.Time`, optional
            Epoch to which to correct refcat proper motion and parallax,
            or `None` to not apply such corrections.

        Returns
        -------
        result : `Photometry` or `Astrometry`
            Result of `fit_function()`
        """
        self.log.info("====== Now processing %s...", name)
        # TODO: this should not print "trying to invert a singular transformation:"
        # if it does that, something's not right about the WCS...
        associations.associateCatalogs(match_cut)
        add_measurement(self.job, 'jointcal.%s_matched_fittedStars' % name,
                        associations.fittedStarListSize())

        applyColorterms = False if name.lower() == "astrometry" else self.config.applyColorTerms
        refCat, fluxField = self._load_reference_catalog(refObjLoader, referenceSelector,
                                                         center, radius, defaultFilter,
                                                         applyColorterms=applyColorterms,
                                                         epoch=epoch)
        refCoordErr = self._get_refcat_coordinate_error_override(refCat, name)

        associations.collectRefStars(refCat,
                                     self.config.matchCut*lsst.geom.arcseconds,
                                     fluxField,
                                     refCoordinateErr=refCoordErr,
                                     rejectBadFluxes=reject_bad_fluxes)
        add_measurement(self.job, 'jointcal.%s_collected_refStars' % name,
                        associations.refStarListSize())

        associations.prepareFittedStars(self.config.minMeasurements)

        self._check_star_lists(associations, name)
        add_measurement(self.job, 'jointcal.%s_prepared_refStars' % name,
                        associations.nFittedStarsWithAssociatedRefStar())
        add_measurement(self.job, 'jointcal.%s_prepared_fittedStars' % name,
                        associations.fittedStarListSize())
        add_measurement(self.job, 'jointcal.%s_prepared_ccdImages' % name,
                        associations.nCcdImagesValidForFit())

        fit_profile_file = 'jointcal_fit_%s.prof'%name if self.config.detailedProfile else ''
        dataName = "{}_{}".format(tract, defaultFilter.physicalLabel)
        with lsst.utils.timer.profile(fit_profile_file):
            result = fit_function(associations, dataName)
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2FilesInitialFinal:
            baseName = self._getDebugPath(f"{name}_final_chi2-{dataName}")
            result.fit.saveChi2Contributions(baseName+"{type}")
            self.log.info("Wrote chi2 contributions files: %s", baseName)

        return result

    def _load_reference_catalog(self, refObjLoader, referenceSelector, center, radius, filterLabel,
                                applyColorterms=False, epoch=None):
        """Load the necessary reference catalog sources, convert fluxes to
        correct units, and apply color term corrections if requested.

        Parameters
        ----------
        refObjLoader : `lsst.meas.algorithms.ReferenceObjectLoader`
            The reference catalog loader to use to get the data.
        referenceSelector : `lsst.meas.algorithms.ReferenceSourceSelectorTask`
            Source selector to apply to loaded reference catalog.
        center : `lsst.geom.SpherePoint`
            The center around which to load sources.
        radius : `lsst.geom.Angle`
            The radius around ``center`` to load sources in.
        filterLabel : `lsst.afw.image.FilterLabel`
            The camera filter to load fluxes for.
        applyColorterms : `bool`
            Apply colorterm corrections to the refcat for ``filterName``?
        epoch : `astropy.time.Time`, optional
            Epoch to which to correct refcat proper motion and parallax,
            or `None` to not apply such corrections.

        Returns
        -------
        refCat : `lsst.afw.table.SimpleCatalog`
            The loaded reference catalog.
        fluxField : `str`
            The name of the reference catalog flux field appropriate for ``filterName``.
        """
        skyCircle = refObjLoader.loadSkyCircle(center,
                                               radius,
                                               filterLabel.bandLabel,
                                               epoch=epoch)

        selected = referenceSelector.run(skyCircle.refCat)
        # Need memory contiguity to get reference filters as a vector.
        if not selected.sourceCat.isContiguous():
            refCat = selected.sourceCat.copy(deep=True)
        else:
            refCat = selected.sourceCat

        if applyColorterms:
            refCatName = refObjLoader.name
            self.log.info("Applying color terms for physical filter=%r reference catalog=%s",
                          filterLabel.physicalLabel, refCatName)
            colorterm = self.config.colorterms.getColorterm(filterLabel.physicalLabel,
                                                            refCatName,
                                                            doRaise=True)

            refMag, refMagErr = colorterm.getCorrectedMagnitudes(refCat)
            refCat[skyCircle.fluxField] = u.Magnitude(refMag, u.ABmag).to_value(u.nJy)
            # TODO: I didn't want to use this, but I'll deal with it in DM-16903
            refCat[skyCircle.fluxField+'Err'] = fluxErrFromABMagErr(refMagErr, refMag) * 1e9

        return refCat, skyCircle.fluxField

    def _check_star_lists(self, associations, name):
        # TODO: these should be len(blah), but we need this properly wrapped first.
        if associations.nCcdImagesValidForFit() == 0:
            raise RuntimeError('No images in the ccdImageList!')
        if associations.fittedStarListSize() == 0:
            raise RuntimeError('No stars in the {} fittedStarList!'.format(name))
        if associations.refStarListSize() == 0:
            raise RuntimeError('No stars in the {} reference star list!'.format(name))

    def _logChi2AndValidate(self, associations, fit, model, chi2Label, writeChi2Name=None):
        """Compute chi2, log it, validate the model, and return chi2.

        Parameters
        ----------
        associations : `lsst.jointcal.Associations`
            The star/reference star associations to fit.
        fit : `lsst.jointcal.FitterBase`
            The fitter to use for minimization.
        model : `lsst.jointcal.Model`
            The model being fit.
        chi2Label : `str`
            Label to describe the chi2 (e.g. "Initialized", "Final").
        writeChi2Name : `str`, optional
            Filename prefix to write the chi2 contributions to.
            Do not supply an extension: an appropriate one will be added.

        Returns
        -------
        chi2: `lsst.jointcal.Chi2Accumulator`
            The chi2 object for the current fitter and model.

        Raises
        ------
        FloatingPointError
            Raised if chi2 is infinite or NaN.
        ValueError
            Raised if the model is not valid.
        """
        if writeChi2Name is not None:
            fullpath = self._getDebugPath(writeChi2Name)
            fit.saveChi2Contributions(fullpath+"{type}")
            self.log.info("Wrote chi2 contributions files: %s", fullpath)

        chi2 = fit.computeChi2()
        self.log.info("%s %s", chi2Label, chi2)
        self._check_stars(associations)
        if not np.isfinite(chi2.chi2):
            raise FloatingPointError(f'{chi2Label} chi2 is invalid: {chi2}')
        if not model.validate(associations.getCcdImageList(), chi2.ndof):
            raise ValueError("Model is not valid: check log messages for warnings.")
        return chi2

    def _fit_photometry(self, associations, dataName=None):
        """
        Fit the photometric data.

        Parameters
        ----------
        associations : `lsst.jointcal.Associations`
            The star/reference star associations to fit.
        dataName : `str`
            Name of the data being processed (e.g. "1234_HSC-Y"), for
            identifying debugging files.

        Returns
        -------
        fit_result : `namedtuple`
            fit : `lsst.jointcal.PhotometryFit`
                The photometric fitter used to perform the fit.
            model : `lsst.jointcal.PhotometryModel`
                The photometric model that was fit.
        """
        self.log.info("=== Starting photometric fitting...")

        # TODO: should use pex.config.RegistryField here (see DM-9195)
        if self.config.photometryModel == "constrainedFlux":
            model = lsst.jointcal.ConstrainedFluxModel(associations.getCcdImageList(),
                                                       self.focalPlaneBBox,
                                                       visitOrder=self.config.photometryVisitOrder,
                                                       errorPedestal=self.config.photometryErrorPedestal)
            # potentially nonlinear problem, so we may need a line search to converge.
            doLineSearch = self.config.allowLineSearch
        elif self.config.photometryModel == "constrainedMagnitude":
            model = lsst.jointcal.ConstrainedMagnitudeModel(associations.getCcdImageList(),
                                                            self.focalPlaneBBox,
                                                            visitOrder=self.config.photometryVisitOrder,
                                                            errorPedestal=self.config.photometryErrorPedestal)
            # potentially nonlinear problem, so we may need a line search to converge.
            doLineSearch = self.config.allowLineSearch
        elif self.config.photometryModel == "simpleFlux":
            model = lsst.jointcal.SimpleFluxModel(associations.getCcdImageList(),
                                                  errorPedestal=self.config.photometryErrorPedestal)
            doLineSearch = False  # purely linear in model parameters, so no line search needed
        elif self.config.photometryModel == "simpleMagnitude":
            model = lsst.jointcal.SimpleMagnitudeModel(associations.getCcdImageList(),
                                                       errorPedestal=self.config.photometryErrorPedestal)
            doLineSearch = False  # purely linear in model parameters, so no line search needed

        fit = lsst.jointcal.PhotometryFit(associations, model)
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2FilesInitialFinal:
            baseName = f"photometry_initial_chi2-{dataName}"
        else:
            baseName = None
        if self.config.writeInitialModel:
            fullpath = self._getDebugPath(f"initial_photometry_model-{dataName}.txt")
            writeModel(model, fullpath, self.log)
        self._logChi2AndValidate(associations, fit, model, "Initialized", writeChi2Name=baseName)

        def getChi2Name(whatToFit):
            if self.config.writeChi2FilesOuterLoop:
                return f"photometry_init-%s_chi2-{dataName}" % whatToFit
            else:
                return None

        # The constrained model needs the visit transform fit first; the chip
        # transform is initialized from the singleFrame PhotoCalib, so it's close.
        if self.config.writeInitMatrix:
            dumpMatrixFile = self._getDebugPath(f"photometry_preinit-{dataName}")
        else:
            dumpMatrixFile = ""
        if self.config.photometryModel.startswith("constrained"):
            # no line search: should be purely (or nearly) linear,
            # and we want a large step size to initialize with.
            fit.minimize("ModelVisit", dumpMatrixFile=dumpMatrixFile)
            self._logChi2AndValidate(associations, fit, model, "Initialize ModelVisit",
                                     writeChi2Name=getChi2Name("ModelVisit"))
            dumpMatrixFile = ""  # so we don't redo the output on the next step

        fit.minimize("Model", doLineSearch=doLineSearch, dumpMatrixFile=dumpMatrixFile)
        self._logChi2AndValidate(associations, fit, model, "Initialize Model",
                                 writeChi2Name=getChi2Name("Model"))

        fit.minimize("Fluxes")  # no line search: always purely linear.
        self._logChi2AndValidate(associations, fit, model, "Initialize Fluxes",
                                 writeChi2Name=getChi2Name("Fluxes"))

        fit.minimize("Model Fluxes", doLineSearch=doLineSearch)
        self._logChi2AndValidate(associations, fit, model, "Initialize ModelFluxes",
                                 writeChi2Name=getChi2Name("ModelFluxes"))

        model.freezeErrorTransform()
        self.log.debug("Photometry error scales are frozen.")

        chi2 = self._iterate_fit(associations,
                                 fit,
                                 self.config.maxPhotometrySteps,
                                 "photometry",
                                 "Model Fluxes",
                                 doRankUpdate=self.config.photometryDoRankUpdate,
                                 doLineSearch=doLineSearch,
                                 dataName=dataName)

        add_measurement(self.job, 'jointcal.photometry_final_chi2', chi2.chi2)
        add_measurement(self.job, 'jointcal.photometry_final_ndof', chi2.ndof)
        return Photometry(fit, model)

    def _fit_astrometry(self, associations, dataName=None):
        """
        Fit the astrometric data.

        Parameters
        ----------
        associations : `lsst.jointcal.Associations`
            The star/reference star associations to fit.
        dataName : `str`
            Name of the data being processed (e.g. "1234_HSC-Y"), for
            identifying debugging files.

        Returns
        -------
        fit_result : `namedtuple`
            fit : `lsst.jointcal.AstrometryFit`
                The astrometric fitter used to perform the fit.
            model : `lsst.jointcal.AstrometryModel`
                The astrometric model that was fit.
            sky_to_tan_projection : `lsst.jointcal.ProjectionHandler`
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

        fit = lsst.jointcal.AstrometryFit(associations, model, self.config.positionErrorPedestal)
        # TODO DM-12446: turn this into a "butler save" somehow.
        # Save reference and measurement chi2 contributions for this data
        if self.config.writeChi2FilesInitialFinal:
            baseName = f"astrometry_initial_chi2-{dataName}"
        else:
            baseName = None
        if self.config.writeInitialModel:
            fullpath = self._getDebugPath(f"initial_astrometry_model-{dataName}.txt")
            writeModel(model, fullpath, self.log)
        self._logChi2AndValidate(associations, fit, model, "Initial", writeChi2Name=baseName)

        def getChi2Name(whatToFit):
            if self.config.writeChi2FilesOuterLoop:
                return f"astrometry_init-%s_chi2-{dataName}" % whatToFit
            else:
                return None

        if self.config.writeInitMatrix:
            dumpMatrixFile = self._getDebugPath(f"astrometry_preinit-{dataName}")
        else:
            dumpMatrixFile = ""
        # The constrained model needs the visit transform fit first; the chip
        # transform is initialized from the detector's cameraGeom, so it's close.
        if self.config.astrometryModel == "constrained":
            fit.minimize("DistortionsVisit", dumpMatrixFile=dumpMatrixFile)
            self._logChi2AndValidate(associations, fit, model, "Initialize DistortionsVisit",
                                     writeChi2Name=getChi2Name("DistortionsVisit"))
            dumpMatrixFile = ""  # so we don't redo the output on the next step

        fit.minimize("Distortions", dumpMatrixFile=dumpMatrixFile)
        self._logChi2AndValidate(associations, fit, model, "Initialize Distortions",
                                 writeChi2Name=getChi2Name("Distortions"))

        fit.minimize("Positions")
        self._logChi2AndValidate(associations, fit, model, "Initialize Positions",
                                 writeChi2Name=getChi2Name("Positions"))

        fit.minimize("Distortions Positions")
        self._logChi2AndValidate(associations, fit, model, "Initialize DistortionsPositions",
                                 writeChi2Name=getChi2Name("DistortionsPositions"))

        chi2 = self._iterate_fit(associations,
                                 fit,
                                 self.config.maxAstrometrySteps,
                                 "astrometry",
                                 "Distortions Positions",
                                 sigmaRelativeTolerance=self.config.astrometryOutlierRelativeTolerance,
                                 doRankUpdate=self.config.astrometryDoRankUpdate,
                                 dataName=dataName)

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
                self.log.warning("ccdImage %s has only %s measuredStars (desired %s)",
                                 ccdImage.getName(), nMeasuredStars, self.config.minMeasuredStarsPerCcd)
            if nRefStars < self.config.minRefStarsPerCcd:
                self.log.warning("ccdImage %s has only %s RefStars (desired %s)",
                                 ccdImage.getName(), nRefStars, self.config.minRefStarsPerCcd)

    def _iterate_fit(self, associations, fitter, max_steps, name, whatToFit,
                     dataName="",
                     sigmaRelativeTolerance=0,
                     doRankUpdate=True,
                     doLineSearch=False):
        """Run fitter.minimize up to max_steps times, returning the final chi2.

        Parameters
        ----------
        associations : `lsst.jointcal.Associations`
            The star/reference star associations to fit.
        fitter : `lsst.jointcal.FitterBase`
            The fitter to use for minimization.
        max_steps : `int`
            Maximum number of steps to run outlier rejection before declaring
            convergence failure.
        name : {'photometry' or 'astrometry'}
            What type of data are we fitting (for logs and debugging files).
        whatToFit : `str`
            Passed to ``fitter.minimize()`` to define the parameters to fit.
        dataName : `str`, optional
            Descriptive name for this dataset (e.g. tract and filter),
            for debugging.
        sigmaRelativeTolerance : `float`, optional
            Convergence tolerance for the fractional change in the chi2 cut
            level for determining outliers. If set to zero, iterations will
            continue until there are no outliers.
        doRankUpdate : `bool`, optional
            Do an Eigen rank update during minimization, or recompute the full
            matrix and gradient?
        doLineSearch : `bool`, optional
            Do a line search for the optimum step during minimization?

        Returns
        -------
        chi2: `lsst.jointcal.Chi2Statistic`
            The final chi2 after the fit converges, or is forced to end.

        Raises
        ------
        FloatingPointError
            Raised if the fitter fails with a non-finite value.
        RuntimeError
            Raised if the fitter fails for some other reason;
            log messages will provide further details.
        """
        if self.config.writeInitMatrix:
            dumpMatrixFile = self._getDebugPath(f"{name}_postinit-{dataName}")
        else:
            dumpMatrixFile = ""
        oldChi2 = lsst.jointcal.Chi2Statistic()
        oldChi2.chi2 = float("inf")
        for i in range(max_steps):
            if self.config.writeChi2FilesOuterLoop:
                writeChi2Name = f"{name}_iterate_{i}_chi2-{dataName}"
            else:
                writeChi2Name = None
            result = fitter.minimize(whatToFit,
                                     self.config.outlierRejectSigma,
                                     sigmaRelativeTolerance=sigmaRelativeTolerance,
                                     doRankUpdate=doRankUpdate,
                                     doLineSearch=doLineSearch,
                                     dumpMatrixFile=dumpMatrixFile)
            dumpMatrixFile = ""  # clear it so we don't write the matrix again.
            chi2 = self._logChi2AndValidate(associations, fitter, fitter.getModel(),
                                            f"Fit iteration {i}", writeChi2Name=writeChi2Name)

            if result == MinimizeResult.Converged:
                if doRankUpdate:
                    self.log.debug("fit has converged - no more outliers - redo minimization "
                                   "one more time in case we have lost accuracy in rank update.")
                    # Redo minimization one more time in case we have lost accuracy in rank update
                    result = fitter.minimize(whatToFit, self.config.outlierRejectSigma,
                                             sigmaRelativeTolerance=sigmaRelativeTolerance)
                    chi2 = self._logChi2AndValidate(associations, fitter, fitter.getModel(), "Fit completed")

                # log a message for a large final chi2, TODO: DM-15247 for something better
                if chi2.chi2/chi2.ndof >= 4.0:
                    self.log.error("Potentially bad fit: High chi-squared/ndof.")

                break
            elif result == MinimizeResult.Chi2Increased:
                self.log.warning("Still some outliers remaining but chi2 increased - retry")
                # Check whether the increase was large enough to cause trouble.
                chi2Ratio = chi2.chi2 / oldChi2.chi2
                if chi2Ratio > 1.5:
                    self.log.warning('Significant chi2 increase by a factor of %.4g / %.4g = %.4g',
                                     chi2.chi2, oldChi2.chi2, chi2Ratio)
                # Based on a variety of HSC jointcal logs (see DM-25779), it
                # appears that chi2 increases more than a factor of ~2 always
                # result in the fit diverging rapidly and ending at chi2 > 1e10.
                # Using 10 as the "failure" threshold gives some room between
                # leaving a warning and bailing early.
                if chi2Ratio > 10:
                    msg = ("Large chi2 increase between steps: fit likely cannot converge."
                           " Try setting one or more of the `writeChi2*` config fields and looking"
                           " at how individual star chi2-values evolve during the fit.")
                    raise RuntimeError(msg)
                oldChi2 = chi2
            elif result == MinimizeResult.NonFinite:
                filename = self._getDebugPath("{}_failure-nonfinite_chi2-{}.csv".format(name, dataName))
                # TODO DM-12446: turn this into a "butler save" somehow.
                fitter.saveChi2Contributions(filename+"{type}")
                msg = "Nonfinite value in chi2 minimization, cannot complete fit. Dumped star tables to: {}"
                raise FloatingPointError(msg.format(filename))
            elif result == MinimizeResult.Failed:
                raise RuntimeError("Chi2 minimization failure, cannot complete fit.")
            else:
                raise RuntimeError("Unxepected return code from minimize().")
        else:
            self.log.error("%s failed to converge after %d steps"%(name, max_steps))

        return chi2

    def _make_output(self, ccdImageList, model, func):
        """Return the internal jointcal models converted to the afw
        structures that will be saved to disk.

        Parameters
        ----------
        ccdImageList : `lsst.jointcal.CcdImageList`
            The list of CcdImages to get the output for.
        model : `lsst.jointcal.AstrometryModel` or `lsst.jointcal.PhotometryModel`
            The internal jointcal model to convert for each `lsst.jointcal.CcdImage`.
        func : `str`
            The name of the function to call on ``model`` to get the converted
            structure. Must accept an `lsst.jointcal.CcdImage`.

        Returns
        -------
        output : `dict` [`tuple`, `lsst.jointcal.AstrometryModel`] or
                 `dict` [`tuple`, `lsst.jointcal.PhotometryModel`]
            The data to be saved, keyed on (visit, detector).
        """
        output = {}
        for ccdImage in ccdImageList:
            ccd = ccdImage.ccdId
            visit = ccdImage.visit
            self.log.debug("%s for visit: %d, ccd: %d", func, visit, ccd)
            output[(visit, ccd)] = getattr(model, func)(ccdImage)
        return output


def make_schema_table():
    """Return an afw SourceTable to use as a base for creating the
    SourceCatalog to insert values from the dataFrame into.

    Returns
    -------
    table : `lsst.afw.table.SourceTable`
        Table with schema and slots to use to make SourceCatalogs.
    """
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    schema.addField("centroid_x", "D")
    schema.addField("centroid_y", "D")
    schema.addField("centroid_xErr", "F")
    schema.addField("centroid_yErr", "F")
    schema.addField("shape_xx", "D")
    schema.addField("shape_yy", "D")
    schema.addField("shape_xy", "D")
    schema.addField("flux_instFlux", "D")
    schema.addField("flux_instFluxErr", "D")
    table = lsst.afw.table.SourceTable.make(schema)
    table.defineCentroid("centroid")
    table.defineShape("shape")
    return table


def get_sourceTable_visit_columns(inColumns, config, sourceSelector):
    """
    Get the sourceTable_visit columns to load from the catalogs.

    Parameters
    ----------
    inColumns : `list`
        List of columns known to be available in the sourceTable_visit.
    config : `JointcalConfig`
        A filled-in config to to help define column names.
    sourceSelector : `lsst.meas.algorithms.BaseSourceSelectorTask`
        A configured source selector to define column names to load.

    Returns
    -------
    columns : `list`
        List of columns to read from sourceTable_visit.
    ixxColumns : `list`
        Name of the ixx/iyy/ixy columns.
    """
    columns = ['visit', 'detector',
               'sourceId', 'x', 'xErr', 'y', 'yErr',
               config.sourceFluxType + '_instFlux', config.sourceFluxType + '_instFluxErr']

    if 'ixx' in inColumns:
        # New columns post-DM-31825
        ixxColumns = ['ixx', 'iyy', 'ixy']
    else:
        # Old columns pre-DM-31825
        ixxColumns = ['Ixx', 'Iyy', 'Ixy']
    columns.extend(ixxColumns)

    if sourceSelector.config.doFlags:
        columns.extend(sourceSelector.config.flags.bad)
    if sourceSelector.config.doUnresolved:
        columns.append(sourceSelector.config.unresolved.name)
    if sourceSelector.config.doIsolated:
        columns.append(sourceSelector.config.isolated.parentName)
        columns.append(sourceSelector.config.isolated.nChildName)
    if sourceSelector.config.doRequireFiniteRaDec:
        columns.append(sourceSelector.config.requireFiniteRaDec.raColName)
        columns.append(sourceSelector.config.requireFiniteRaDec.decColName)
    if sourceSelector.config.doRequirePrimary:
        columns.append(sourceSelector.config.requirePrimary.primaryColName)

    return columns, ixxColumns


def extract_detector_catalog_from_visit_catalog(table, visitCatalog, detectorId,
                                                ixxColumns, sourceFluxType, log):
    """Return an afw SourceCatalog extracted from a visit-level dataframe,
    limited to just one detector.

    Parameters
    ----------
    table : `lsst.afw.table.SourceTable`
        Table factory to use to make the SourceCatalog that will be
        populated with data from ``visitCatalog``.
    visitCatalog : `pandas.DataFrame`
        DataFrame to extract a detector catalog from.
    detectorId : `int`
        Numeric id of the detector to extract from ``visitCatalog``.
    ixxColumns : `list` [`str`]
        Names of the ixx/iyy/ixy columns in the catalog.
    sourceFluxType : `str`
        Name of the catalog field to load instFluxes from.
    log : `logging.Logger`
        Logging instance to log to.

    Returns
    -------
    catalog : `lsst.afw.table.SourceCatalog`, or `None`
        Detector-level catalog extracted from ``visitCatalog``, or `None`
        if there was no data to load.
    """
    # map from dataFrame column to afw table column
    mapping = {'x': 'centroid_x',
               'y': 'centroid_y',
               'xErr': 'centroid_xErr',
               'yErr': 'centroid_yErr',
               ixxColumns[0]: 'shape_xx',
               ixxColumns[1]: 'shape_yy',
               ixxColumns[2]: 'shape_xy',
               f'{sourceFluxType}_instFlux': 'flux_instFlux',
               f'{sourceFluxType}_instFluxErr': 'flux_instFluxErr',
               }

    catalog = lsst.afw.table.SourceCatalog(table)
    matched = visitCatalog['detector'] == detectorId
    n = sum(matched)
    if n == 0:
        return None
    catalog.resize(sum(matched))
    view = visitCatalog.loc[matched]
    catalog['id'] = view.index
    for dfCol, afwCol in mapping.items():
        catalog[afwCol] = view[dfCol]

    log.debug("%d sources selected in visit %d detector %d",
              len(catalog),
              view['visit'].iloc[0],  # all visits in this catalog are the same, so take the first
              detectorId)
    return catalog
