# Render recipes

The scene files and render scripts for the published gallery. The rendered
output itself is not in this repository: it lives at

  https://dreevich.dev/gr_raytracer-gallery/   (mdreem/gr_raytracer-gallery)

with the full-resolution PNGs, the linear `.hdr` masters behind each grade, and
the animations in object storage, linked under every figure there. The captions
and the per-frame grade parameters that used to sit in `images/*.md` moved to
that repository too.

What stays here:

- `*.toml`, the scene definitions each gallery frame was rendered from
- `create-*.sh`, the scripts that drive batches of them
- `sky-survey/sky-survey.toml`, the flat-space scene used to survey the sky
- the handful of PNGs the top-level `README.md` embeds

Rendering a frame from one of these scenes writes a linear `.hdr`; grading it to
a PNG is a separate step (`scripts/grade.py`). The gallery records which white
point, bloom and tone curve each published frame used.
