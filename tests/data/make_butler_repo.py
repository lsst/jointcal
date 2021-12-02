#!/bin/env python
# Make a gen3 butler and import the files to have a working cfht_minimal repo.
#
# This can be used to regenerate the gen3 registry when the butler has changes
# that are not backwards compatible.
#
# I had to temporarily alter the sdss-dr9-fink-v5b.ecsv file to have absolute
# paths, because I couldn't get `ingest_files` to work with `prefix=` no matter
# what I did on this script, even though it works for the tests.

import lsst.daf.butler
from lsst.daf.butler.script import ingest_files
import lsst.obs.base

repopath = "cfht_minimal/repo"
instrument = "lsst.obs.cfht.MegaPrime"
lsst.daf.butler.Butler.makeRepo(repopath)
butler = lsst.daf.butler.Butler(repopath, writeable=True)
instrInstance = lsst.obs.base.utils.getInstrument(instrument)
instrInstance.register(butler.registry)

graph = butler.registry.dimensions.extract(["htm7"])
datasetType = lsst.daf.butler.DatasetType("sdss_dr9_fink_v5b", graph,
                                          "SimpleCatalog", universe=butler.registry.dimensions)
butler.registry.registerDatasetType(datasetType)
ingest_files(repopath, "sdss_dr9_fink_v5b", "refcats",
             "cfht_minimal/sdss-dr9-fink-v5b.ecsv")
butler.import_(directory=repopath, filename="cfht_minimal/exports.yaml",
               skip_dimensions={'instrument', 'detector', 'physical_filter'})
