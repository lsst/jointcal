# jointcal/tests/data

## cfht_minimal

A trivial dataset to do a minimal test of jointcal functionality, without having to download the full testdata_jointcal repo.
This data is small enough that the final photometry chi2 can be computed a-priori by solving the matrix system directly.

Useful for debugging photometry code. A "minimal" catalog extracted from cfht (see above), containing 2 ccds with 2 sources each.
Two sources in one catalog have a refcat match (using the `cfht` `sdss-dr9-fink-v5b` reference catalog), and one source is matched between the catalogs.
Contains `ccd=12` with (0-indexed) `rows=336,337` of `visit=850587` and `rows=139,140` for `visit=849375`, all taken from the `cfht/` repo in [testdata_jointcal](https://github.com/lsst/testdata_jointcal)
In jointcal, the photometry fit will contain 3 valid measuredStars, 2 fittedStars, and 2 refStars.

## parquet test data

``extractParquetSubset.py`` trims a ``sourceTable_visit`` parquet catalog from ``testdata_jointcal/hsc`` to use to test jointcal's code to covert a visit-level dataframe into multiple detector-level afw tables.
