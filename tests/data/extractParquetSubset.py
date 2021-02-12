#!/usr/bin/env python
"""Cut down a parquet source table to use in unittests.
"""
import os

import numpy as np
import pyarrow.parquet

filename = f"{os.environ['TESTDATA_JOINTCAL_DIR']}/hsc/01291/HSC-R/output/sourceTable-0034690.parq"
file = pyarrow.parquet.ParquetFile(filename)
data = file.read(use_pandas_metadata=True).to_pandas()

ccds = set(data['ccd'])
selection = []
for ccd in ccds:
    # extract a number of elements from each detector-catalog equal to that ccd id
    subsample = np.where(data['ccd'] == ccd)[0][:ccd]
    selection.extend(subsample)

subset = data.iloc[selection]
parq = pyarrow.Table.from_pandas(subset)
outfile = "subselected-sourceTable-0034690.parq"
pyarrow.parquet.write_table(parq, outfile, compression="none")
n = 100
print(f"Row {n}; ccd, id, x, y")
print(subset.ccd.iloc[n], subset.sourceId.iloc[n], subset.x.iloc[n], subset.y.iloc[n])
print(f"Wrote ccds {ccds} into: {outfile}")
