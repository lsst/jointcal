from lsst.pipe.tasks.colorterms import ColortermLibrary
# This file is to test that Config validation works when colorterms are
# requested but no ColortermLibrary is provided.
config.applyColorTerms = True
# Have to override the fact that the default HSC jointcal config provides
# a ColortermLibrary.
config.colorterms = ColortermLibrary()
