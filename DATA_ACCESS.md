# Data Access

The empirical application uses licensed 1-minute S&P 500 data from TickData.

The raw TickData file is not included in this repository and should not be committed. To run the empirical scripts,
place the licensed file at:

```text
empirical_analysis/data/SPX.csv
```

or pass the file explicitly:

```sh
Rscript empirical_analysis/run_empirical_bundle_refresh.R --data_file=/path/to/SPX.csv
```

The expected CSV contains at least:

- `DateTime`: timestamp string readable as `%Y-%m-%d %H:%M:%S`
- `Close`: intraday S&P 500 price level used to construct one-minute log-return curves

Access to the raw data is governed by the relevant TickData license. Users who do not have licensed access can still
run the critical-value and simulation code.
