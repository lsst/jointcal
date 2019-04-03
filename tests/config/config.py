config.astrometryRefObjLoader.ref_dataset_name = "sdss-dr9-fink-v5b"
config.photometryRefObjLoader.ref_dataset_name = "sdss-dr9-fink-v5b"
# TODO DM-17597: Use the astrometrySourceSelector until we redo the refcats and metrics.
# Once we do that, we can try using the sourceSelector defaults and recompute all
# the metrics to match, though it may result in problems due to not having enough sources
# because the testdata is generally fainter than the HSC/LSST data that the defaults
# are now designed for.
config.sourceSelector.name = "astrometry"
config.sourceSelector["astrometry"].sourceFluxType = "Calib"
config.sourceSelector["astrometry"].badFlags.extend(["slot_Shape_flag", "base_PixelFlags_flag_interpolated"])
config.astrometryReferenceErr = 100  # milli-arcsecond
