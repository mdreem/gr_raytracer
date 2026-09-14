#!/usr/bin/env -S uv run --group gaia python
"""Download a magnitude-limited subset of the Gaia DR3 star catalogue.

Queries ESA's Gaia TAP/ADQL service server-side and saves only the requested
subset locally as a typed Parquet file. The renderer loads this file to
compute celestial directions, fluxes, and colour temperatures. Rendering
and spatial indexing are outside this script's scope.

Usage:
    uv run --group gaia scripts/gaia/download.py download --max-magnitude 12 --limit 100000
    uv run --group gaia scripts/gaia/download.py download --max-magnitude 12

Run from the repository root.
"""

from __future__ import annotations

import argparse
import os
import sys

# Columns retrieved from gaiadr3.gaia_source, in query/output order.
GAIA_TABLE = "gaiadr3.gaia_source"
COLUMNS = [
    "source_id",
    "ra",
    "dec",
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
    "bp_rp",
]

# Photometry the renderer needs; rows missing any of these are useless as
# point sources, so we require them to be non-null server-side.
REQUIRED_NON_NULL = [
    "phot_g_mean_mag",
    "phot_bp_mean_mag",
    "phot_rp_mean_mag",
]

# Parquet dtypes. Gaia source_id is a 64-bit designation that fits in int64
# (max ~6.9e18 < 9.2e18); positions want full float64 precision; photometry
# and colour are fine as float32 to keep the file compact.
DTYPES = {
    "source_id": "int64",
    "ra": "float64",
    "dec": "float64",
    "phot_g_mean_mag": "float32",
    "phot_bp_mean_mag": "float32",
    "phot_rp_mean_mag": "float32",
    "bp_rp": "float32",
}

DEFAULT_OUTPUT = "data/gaia_dr3.parquet"


def build_adql_query(max_magnitude: float, limit: int | None) -> str:
    """Build the ADQL query string sent to the Gaia TAP service.

    Only a ``TOP`` clause is added when ``limit`` is given; an omitted limit
    means the full magnitude-limited result set is returned.
    """
    top = f"TOP {limit} " if limit is not None else ""
    select_cols = ", ".join(COLUMNS)
    non_null = "\n  AND ".join(f"{c} IS NOT NULL" for c in REQUIRED_NON_NULL)
    return (
        f"SELECT {top}{select_cols}\n"
        f"FROM {GAIA_TABLE}\n"
        f"WHERE phot_g_mean_mag <= {max_magnitude}\n"
        f"  AND {non_null}"
    )


def _human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024.0 or unit == "TiB":
            return f"{size:.1f} {unit}"
        size /= 1024.0
    return f"{size:.1f} TiB"


def download(max_magnitude: float, limit: int | None, output: str, overwrite: bool) -> None:
    # Imported lazily so `--help` works without the astro stack installed.
    import pyarrow as pa
    import pyarrow.parquet as pq
    from astroquery.gaia import Gaia

    if os.path.exists(output) and not overwrite:
        sys.exit(
            f"Output {output!r} already exists. Use --overwrite to replace it."
        )

    query = build_adql_query(max_magnitude, limit)
    print("ADQL query:")
    print(query)
    print()

    # flush=True so progress is visible during the long async job even when
    # stdout is a pipe (Python block-buffers non-tty output otherwise).
    print("Launching asynchronous Gaia job ...", flush=True)
    job = Gaia.launch_job_async(query)
    table = job.get_results()

    # astropy Table -> pandas; masked/absent photometry becomes NaN. The
    # server-side filter already drops null G/BP/RP, but guard locally too.
    df = table.to_pandas()
    df = df.dropna(subset=REQUIRED_NON_NULL)

    # Coerce to the compact, precision-preserving dtypes.
    df = df[COLUMNS].astype(DTYPES)

    print(f"Returned {len(df)} stars.", flush=True)

    out_dir = os.path.dirname(output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    arrow_table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(arrow_table, output, compression="zstd")

    size = os.path.getsize(output)
    print(f"Wrote {output} ({_human_size(size)}, {size} bytes).")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Gaia DR3 star-catalogue helper for the raytracer.",
    )
    # Subcommands so a later `preprocess` stage can slot in alongside `download`.
    sub = parser.add_subparsers(dest="command", required=True)

    dl = sub.add_parser(
        "download",
        help="Download a magnitude-limited Gaia DR3 subset to Parquet.",
    )
    dl.add_argument(
        "--max-magnitude",
        type=float,
        default=12.0,
        help="Keep stars with phot_g_mean_mag <= this value (default: 12).",
    )
    dl.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on returned rows (adds ADQL TOP); handy for testing.",
    )
    dl.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Output Parquet path (default: {DEFAULT_OUTPUT}).",
    )
    dl.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "download":
        download(
            max_magnitude=args.max_magnitude,
            limit=args.limit,
            output=args.output,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    main()
