#!/usr/bin/env python

"""
Script used to rename columns in the local test source table parquet files,
in particular renaming gen2 ``ccd`` to modern ``detector`` from DM-31889.
"""

import os
import glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import lsst.utils


data_dir = os.path.join(lsst.utils.getPackageDir("jointcal"), "tests", "data")

directories = [
    "cfht_minimal/repo/src/06AL01/D3/2006-05-20/r",
    "cfht_minimal/repo/src/06AL01/D3/2006-06-02/r",
    "",
]

for directory in directories:
    parquet_files = sorted(glob.glob(os.path.join(data_dir, directory, "*.parq")))

    for parquet_file in parquet_files:
        df = pd.read_parquet(parquet_file)

        columns = list(df.columns)

        mapper = {}
        for column in columns:
            if column == "ccd":
                new_column = "detector"
                mapper[column] = new_column

        df = df.rename(columns=mapper)

        os.remove(parquet_file)

        table = pa.Table.from_pandas(df)
        pq.write_table(table, parquet_file)
