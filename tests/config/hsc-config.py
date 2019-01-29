# no colorterms in tests unless specifically requested (this overrides the default HSC config)
config.applyColorTerms = False
config.photometryRefObjLoader.filterMap = {'y': 'z'}
config.astrometryRefObjLoader.filterMap = {'y': 'z'}
