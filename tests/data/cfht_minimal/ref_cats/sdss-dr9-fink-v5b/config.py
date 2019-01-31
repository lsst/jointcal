import lsst.meas.algorithms.ingestIndexReferenceTask
assert type(config)==lsst.meas.algorithms.ingestIndexReferenceTask.DatasetConfig, 'config is of type %s.%s instead of lsst.meas.algorithms.ingestIndexReferenceTask.DatasetConfig' % (type(config).__module__, type(config).__name__)
# String to pass to the butler to retrieve persisted files.
config.ref_dataset_name='sdss'

# Depth of the HTM tree to make.  Default is depth=7 which gives
#               ~ 0.3 sq. deg. per trixel.
config.indexer['HTM'].depth=7

config.indexer.name='HTM'