# Joint Calibration

This package performs joint astrometry and photometry across multiple visits of the same field, fitting models for the astrometric distortions and photometric calibrations of the CCDs and focal plane.

It was derived from lsst_france/meas_simastrom

## A note on tests

All of the `tests/test_jointcal_<obs>.py` tests run jointcal on their respective `testdata_jointcal/<obs>` dataset.
The numeric values being tested (the metrics, pa1, and relative/absolute rms) were not computed in advance, but are empirical, from the first time that test was successfully run.
They are updated when changes to the algorithm warrants, e.g. a change to outlier rejection might affect the final chi2/ndof for some of the testdata.
More sophisticated testing of jointcal's validity is handled by `validate_drp` on much larger datasets that have multiple full-focal plane visits.
