# Joint Calibration

This package performs astrometry and photometry across multiple visits of the same field, fitting models for the astrometric distortions and photometric calibrations of the CCDs and focal plane.

It was derived from [meas_simastrom](https://github.com/lsst-france/meas_simastrom).

## A note on numbers in tests

All of the `tests/test_jointcal_<obs>.py` tests run jointcal on their respective `testdata_jointcal/<obs>` dataset.
The numeric values being tested in each test's `metrics` dictionary were not computed in advance, but are empirical, from the first time that test was successfully run.
The one exception is the `cfht_minimal` dataset, which was designed to be small enough (3 `MeasuredStar`, 2 `FittedStar`, 1 `RefStar`) that the "correct" answers could be computed in advance with e.g. Mathematica.
They are updated when changes to the algorithm warrants, e.g. a change to outlier rejection that affects the final chi2/ndof for some of the testdata.

More sophisticated testing of jointcal's validity is handled by `faro` on much larger datasets that have multiple full-focal plane visits.
