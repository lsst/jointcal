# jointcal/tests/data

## cfht_minimal

A trivial dataset to do a minimal test of jointcal functionality, without having to download the full testdata_jointcal repo.
This data is small enough that the final photometry chi2 can be computed a-priori by solving the matrix system directly.

Useful for debugging photometry code. A "minimal" catalog extracted from cfht (see above), containing 2 ccds with 2 sources each, with the visits centered at `214.884832, 52.6622199`.
Two sources in one catalog have a refcat match (using the `cfht` `sdss-dr9-fink-v5b` reference catalog), and one source is matched between the catalogs.
Contains `ccd=12` with (0-indexed) `rows=336,337` of `visit=850587` and `rows=139,140` for `visit=849375`, all taken from the `cfht/` repo in [testdata_jointcal](https://github.com/lsst/testdata_jointcal)
In jointcal, the photometry fit will contain 3 valid measuredStars, 2 fittedStars, and 2 refStars.

It also contains a working gen3 butler repo, generated with `make_butler_repo.py`, which allows `testUtils.py` to use this data to construct simulated ccdImages.
The `exports.yaml` file (used to create that gen3 repo and to create temporary repos in `test_jointcal_chft_minimal.py`) was hand-edited from the one included in the `testdata_jointcal` package's `cfht/` directory; see that package for how that file was created.

The parquet catalogs were updated to replace gen2 ``ccd`` with gen3 ``detector`` on 25-Oct-2022 with the script `rename_sourcetable_columns.py`.

## parquet test data

``extractParquetSubset.py`` trims a ``sourceTable_visit`` parquet catalog from ``testdata_jointcal/hsc`` to use to test jointcal's code to covert a visit-level dataframe into multiple detector-level afw tables.

## output

``HSC-R/9615`` contains a subset of the output from the ``RC/w_2021_10/DM-29074`` gen2 run on lsst-devl, containing 3 visits (23900, 23924, 28976) with four detectors each (37, 49, 90, 103).
This data is used to test the cameraGeometry calculation code in jointcal.
