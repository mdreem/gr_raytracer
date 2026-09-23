# Generated data

This directory holds **generated** catalogue data used by the renderer. It is
tracked (via this README) only so the directory always exists in a checkout,
so a build or run can rely on `data/` being present. Everything else in here is
git-ignored and must not be committed (see the repo `.gitignore`).

## Gaia DR3 star catalogue

The star-catalog helper (`docs/plan-12-star-catalog.md`) writes its Parquet
download here by default:

```sh
uv run --group gaia scripts/gaia/download.py download --max-magnitude 12 --output data/gaia_dr3.parquet
```

See `scripts/gaia/README.md` for options and the output schema. The renderer
loads this Parquet file directly through the `[star_catalog]` scene config;
there is no separate binary-conversion step.
