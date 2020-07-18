# no colorterms in tests unless specifically requested (this overrides the default HSC config)
config.applyColorTerms = False

# Use Gaia-DR2 for astrometry (override HSC jointcal defaults, which are PS1 for both)
config.astrometryRefObjLoader.ref_dataset_name = "gaia_dr2_20200414"
config.astrometryRefObjLoader.filterMap = {}
config.astrometryRefObjLoader.anyFilterMapsToThis = 'phot_g_mean'
config.astrometryReferenceErr = None
