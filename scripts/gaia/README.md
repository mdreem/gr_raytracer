# Gaia DR3 star-catalogue helper

Query ESA's Gaia TAP/ADQL service server-side and save a magnitude-limited
subset locally as a typed Parquet file (see `docs/plan-12-star-catalog.md`).
Only the requested subset is downloaded, not the Gaia bulk archive.

The renderer loads this Parquet file directly: point it at the file with the
`path` key of the `[star_catalog]` section in a scene TOML. There is no
separate binary-conversion step.

## Dependencies

The astro stack (`astroquery`, `astropy`, `pyarrow`) lives in an opt-in
`gaia` dependency group in the root `pyproject.toml`, so it is not pulled
in by a plain `uv sync`. `uv run --group gaia ...` installs it on first use.

## Usage

Run from the repository root.

Small test download (~100k stars, `TOP` clause added):

```sh
uv run --group gaia scripts/gaia/download.py download --max-magnitude 12 --limit 100000
```

Larger, full magnitude-limited query (no `TOP`, all stars brighter than G=12):

```sh
uv run --group gaia scripts/gaia/download.py download --max-magnitude 12
```

### Options

| Option            | Default                  | Meaning                                        |
| ----------------- | ------------------------ | ---------------------------------------------- |
| `--max-magnitude` | `12`                     | Keep stars with `phot_g_mean_mag <=` this.     |
| `--limit`         | (none)                   | Cap rows via ADQL `TOP`; omit for the full set.|
| `--output`        | `data/gaia_dr3.parquet`  | Output Parquet path.                           |
| `--overwrite`     | off                      | Overwrite the output file if it exists.        |

## Output format

Parquet (Zstandard compression) with one row per star:

| Column             | Type    | Notes                              |
| ------------------ | ------- | ---------------------------------- |
| `source_id`        | int64   | Gaia DR3 designation.              |
| `ra`, `dec`        | float64 | Degrees (ICRS).                    |
| `phot_g_mean_mag`  | float32 | G-band mean magnitude.             |
| `phot_bp_mean_mag` | float32 | BP-band mean magnitude.            |
| `phot_rp_mean_mag` | float32 | RP-band mean magnitude.            |
| `bp_rp`            | float32 | BP - RP colour.                    |

Rows missing G, BP, or RP photometry are filtered out server-side (and
again locally as a guard). Generated catalogue files under `data/` are
git-ignored and should not be committed.

## Read it back

```python
import pandas as pd
df = pd.read_parquet("data/gaia_dr3.parquet")
print(df.dtypes, len(df))
```
