# To test that jointcal fails when loading a refcat without coordinate errors
# and with no "artificial error" suplied.
# We have to use the sdss refcat, as gaia supplies coord errors.
config.astrometryRefObjLoader.ref_dataset_name = "sdss-dr9-fink-v5b"
config.astrometryRefObjLoader.filterMap = {}
config.astrometryReferenceErr = None
