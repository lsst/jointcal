config.sourceSelector.name = "flagged"
# Calib flag names changed with RFC-498 (DM-14997).  The following sets the config to use the
# old names associated with the current data in testdata_jointcal that was processed pre-RFC-498.
# Remove line if the data in testdata_jointcal are ever reprocessed post-RFC-498 (e.g. DM-17597)
config.sourceSelector.active.field = "calib_psfUsed"
