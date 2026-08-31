# Generated data

This directory holds **generated** catalogue data used by the renderer. It is
tracked (via this README) only so the directory always exists in a checkout,
so a build or run can rely on `data/` being present. Everything else in here is
git-ignored and must not be committed (see the repo `.gitignore`).

## Gaia DR3 star catalogue

Stage 1 of the star-catalog pipeline (`docs/plan-12-star-catalog.md`) writes
its Parquet download here by default:

```sh
uv run --group gaia scripts/gaia/download.py download --max-magnitude 12 --output data/gaia_dr3.parquet
```

See `scripts/gaia/README.md` for options and the output schema. A later
preprocessing stage will convert this Parquet into the renderer's compact
binary point-source format, also written under `data/`.
