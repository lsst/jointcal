#!/usr/bin/bash
#
# Process the most recent validation_data_hsc dataset.
# Requires the following be setup before being run:
#     obs_subaru
#     jointcal
#     validation_data_hsc

OUTPUT="hsc_output"

VISITS="903982^904006^904828^904846"
ID="--id visit=$VISITS"

CLOBBER="--clobber-versions --clobber-config"

PROFILE="--profile jointcal-hsc.prof"

echo jointcal.py $VALIDATION_DATA_HSC_DIR/DATA/rerun/20160805 --output $OUTPUT $ID $CLOBBER $PROFILE
jointcal.py $VALIDATION_DATA_HSC_DIR/DATA/rerun/20160805 --output $OUTPUT $ID $CLOBBER $PROFILE
