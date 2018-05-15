# Joint Calibration

This package performs joint astrometry and photometry across multiple visits of the same field, fitting models for the astrometric distortions and photometric calibrations of the CCDs and focal plane.

It was derived from lsst_france/meas_simastrom

## A note on numbers in tests

All of the `tests/test_jointcal_<obs>.py` tests run jointcal on their respective `testdata_jointcal/<obs>` dataset.
The numeric values being tested (the `metrics` dictionary, and the `pa1`/`relative_error`/`dist_absolute_rms` values) were not computed in advance, but are empirical, from the first time that test was successfully run.
The one exception is the `cfht_minimal` dataset, which was designed to be small enough (3 `MeasuredStar`, 2 `FittedStar`, 1 `RefStar`) that the "correct" answers could be computed in advance with e.g. Mathematica.
They are updated when changes to the algorithm warrants, e.g. a change to outlier rejection that affects the final chi2/ndof for some of the testdata.
The `pa1` (see LPM-17 table 14) and relative/absolute rms (related to, but not the same as LPM-17's AM1 in table 18) values should be considered regression tests: algorithmic changes shouldn't result in larger values, but could produce smaller ones.
They are "measured value must be less than" tests, and the calculations are in `lsst.jointcal.utils.JointcalStatistics.py`.

More sophisticated testing of jointcal's validity is handled by `validate_drp` on much larger datasets that have multiple full-focal plane visits.
