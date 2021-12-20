#!/bin/env python
# Generate the visit and source table summary files for the cfht_minimal
# dataset.
# This was run with a hand-modified source.yaml file, with things removed
# to allow this quite old data to be converted to minimally-useful parqet files.
# This was written using the gen2 postprocess `runDataRef` calls, as it was
# part of the process of converting this repo to gen3.

from lsst.pipe.tasks import postprocess
import lsst.daf.persistence

import logging
logging.basicConfig(level=logging.DEBUG, format="{name} {levelname}: {message}", style="{")

butler = lsst.daf.persistence.Butler('cfht_minimal/repo')

# each visit has to be run separately
ref1 = butler.dataRef('calexp', visit=849375, ccd=12)
ref2 = butler.dataRef('calexp', visit=850587, ccd=12)

task = postprocess.ConsolidateVisitSummaryTask()
task.runDataRef([ref1])
task.runDataRef([ref2])

task = postprocess.WriteSourceTableTask()
task.runDataRef(ref1)
task.runDataRef(ref2)

# TransformSourceTableTask needs the parquet source files as dataRefs.
ref1 = butler.dataRef('source', visit=849375, ccd=12)
ref2 = butler.dataRef('source', visit=850587, ccd=12)
config = postprocess.TransformSourceTableTask.ConfigClass()
config.functorFile = "cfht_minimal_sourceTable.yaml"
task = postprocess.TransformSourceTableTask(config=config)
task.runDataRef(ref1)
task.runDataRef(ref2)

# ConsolidateSourceTableTask needs the parquet sourceTable as dataRefs.
ref1 = butler.dataRef('sourceTable', visit=849375, ccd=12)
ref2 = butler.dataRef('sourceTable', visit=850587, ccd=12)
task = postprocess.ConsolidateSourceTableTask()
task.runDataRef([ref1])
task.runDataRef([ref2])
