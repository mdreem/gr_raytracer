# Starfield gallery recipes

The blackbody, volumetric, and temperature-series frames in
[`kerr.md`](kerr.md) were re-rendered over the **real Gaia DR3 star
catalogue** (point sources to magnitude 12, `data/gaia_mag12.parquet`),
replacing the earlier image-texture backgrounds (M25, NGC 6355). Each scene
adds a `[star_catalog]` block so the catalogue stars are lensed and
redshifted along with the rays, then the linear HDR is graded to PNG with
[`scripts/grade.py`](../scripts/grade.py).

Two-step pipeline per frame:

1. **Render** the scene to a linear `.hdr`. Common flags:
   `--max-steps=1000000 --epsilon=1e-9 --theta=-3.14159 --psi=0 --phi=0`
   (offloaded via `scripts/render-remote.sh` for the heavier frames).
2. **Grade**: `uv run scripts/grade.py in.hdr out.png --white W --bloom B --tonemap {aces,reinhard}`.

The grade white point sits at roughly the disc's 99th-percentile luminance
divided by ~227, which keeps the disc contained (near white) while the much
fainter catalogue stars land in visible mid-tones.

## Blackbody discs (Novikov-Thorne, `scene-definitions/kerr-blackbody-disc.toml` + `[star_catalog]`)

Scene: `flux_scale = 8000`, `celestial_temperature = 0` (catalogue stars only).

| image | a | T (K) | inner_radius | camera | res | white | bloom | tonemap |
|-------|---|-------|--------------|--------|-----|-------|-------|---------|
| `kerr_blackbody_disk_1` (gold far) | 0.499 | 8000 | 0.63 | -22,0,0.9 | 1000² | 450 | 0.08 | aces |
| `kerr_blackbody_disk_2` (gold near) | 0.499 | 8000 | 0.63 | -6,0,0.45 | 1000² | 50000 | 0.08 | aces |
| `kerr_blackbody_disk_a_0_5__1` (red far) | 0.4995 | 4000 | 0.63 | -22,0,0.9 | 1000² | 0.2 | 0.08 | aces |
| `kerr_blackbody_disk_a_0_5__2` (red near, high-contrast) | 0.4995 | 4000 | 0.603 | -6,0,0.45 | 1000² | 1000 | 0.08 | aces |
| `kerr_blackbody_disk_a_0_5__2_lifted` (red near, lifted) | 0.4995 | 4000 | 0.603 | -6,0,0.45 | 1000² | 200 | 0.05 | reinhard |

`inner_radius` tracks `1.02 * r_isco(a)`: 0.63 at a = 0.499, 0.603 at the
essentially-extremal a = 0.4995. The lifted red-near frame is the same HDR as
the high-contrast one, graded with a gentler Reinhard toe so the cool outer
disc stays visible.

## Temperature series (`brightness_reference_temperature = 7000`, inner 0.8, a = 0.499)

Scene: the main-image volumetric disc + `[star_catalog]`, camera `-17,0,1.5`,
1280x720. Only the peak temperature and the brightness normalisation change.

| image | T (K) | flux_scale | white | bloom | tonemap |
|-------|-------|-----------|-------|-------|---------|
| `kerr_disc_temperature_8000k` | 8000 | 20 | 13 | 0.1 | aces |
| `kerr_disc_temperature_12000k` | 12000 | 1000 | 735 | 0.1 | aces |
| `kerr_disc_temperature_20000k` | 20000 | 70000 | 47000 | 0.1 | aces |

## Volumetric (`scene-definitions/kerr-bl-volumetric-streaky.toml` + `[star_catalog]`)

| image | camera | res | white | bloom | tonemap |
|-------|--------|-----|-------|-------|---------|
| `render_kerr_stars_volumetric` | -20,0,-0.6 | 1500² | 0.5 | 0.12 | aces |
