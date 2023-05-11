# Config overrides for jointcal tests.

# jointcal test metrics were made including unresolved sources.
config.sourceSelector['science'].doUnresolved = False
# Test catalogs do not have detect_isPrimary
config.sourceSelector['science'].doRequirePrimary = False
