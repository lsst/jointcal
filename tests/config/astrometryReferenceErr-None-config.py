# To test that jointcal fails when loading a refcat without coordinate errors
# and with no "artificial error" suplied.
config.astrometryReferenceErr = None

# We have to use the sdss refcat, as gaia includes coord errors; this tests
# how jointcal handles the case where there is no coord error in the  refcat
# and no `astrometryReferenceErr` is supplied (should get a useful error message).
config.connections.astrometryRefCat = "sdss_dr9_fink_v5b"
# jointcal overrides the refcat loader default for gaia, so we have to clear it for SDSS.
config.astrometryRefObjLoader.anyFilterMapsToThis = None
# This test refcat cannot apply proper motions.
config.astrometryRefObjLoader.requireProperMotion = False
# This old SDSS catalog has capitalized filter names, for some reason.
config.astrometryRefObjLoader.filterMap = {'r': 'R'}
