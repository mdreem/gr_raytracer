# Schwarzschild Black Hole

[← Back to the gallery index](images.md)

### Gaia star field

A Schwarzschild black hole lensing the real Gaia DR3 catalogue (G &le; 12) as
point-source stars, no accretion disc. The shadow and photon ring are perfectly
circular; the ring glow is clumpy only because the real sky is uneven (recipe:
`images/gaia-mag12-schwarzschild-no-disc-1920x1080-2026-09-21.toml`; rendered to
linear HDR and graded with bloom + a luminance-preserving ACES tone map at
exposure 86).

<div align="center">
  <img src="gaia-mag12-schwarzschild-no-disc-1920x1080-2026-09-21.png" alt="A Schwarzschild black hole gravitationally lensing the real Gaia DR3 star field">
  <p>Gaia DR3 stars (G &le; 12) lensed by a Schwarzschild black hole, viewed edge-on from the &minus;x axis; 1920&times;1080 HDR, bloom + ACES.</p>
</div>

The same star field from the vantage the disc images use (camera `0,-24,6.5`,
looking back at the hole), so the lensed sky patch matches those below (recipe:
`images/gaia-mag12-schwarzschild-no-disc-discview-1920x1080-2026-09-23.toml`;
same grade and exposure).

<div align="center">
  <img src="gaia-mag12-schwarzschild-no-disc-discview-1920x1080-2026-09-23.png" alt="The Gaia DR3 star field lensed by a Schwarzschild black hole, seen from the accretion-disc vantage point">
  <p>Same catalogue from the disc images' vantage: the sky patch here is the one the disc renders below sit in front of.</p>
</div>

### Accretion disc over the Gaia star field

A flat Novikov-Thorne blackbody accretion disc rendered together with the real
Gaia DR3 (G &le; 12) star field in a single pass. The disc's `flux_scale` dims
its emitted light into the star field's brightness range while leaving its
opacity untouched, so the disc still fully occludes the stars behind it and one
exposure captures both. Rendered to linear HDR and graded with bloom and a
luminance-preserving ACES tone map (grading brightness while keeping
chromaticity, so the disc holds its colour instead of washing to white). Same
camera and star field in both; only the disc temperature differs
(recipes: `images/schwarzschild-disc-stars-yellow-1920x1080-2026-09-23.toml`,
`images/schwarzschild-disc-stars-blue-1920x1080-2026-09-23.toml`; exposure 86).
An HDR gain-map JPEG sits beside each PNG (`*-hdr.jpg`): an ordinary SDR JPEG
carrying a hidden gain map, so it shows the SDR image everywhere and brightens
the disc past white on an HDR screen (built with the macOS CoreImage HDR API).

<div align="center">
  <img src="schwarzschild-disc-stars-yellow-1920x1080-2026-09-23.png" alt="A warm 4500 K blackbody accretion disc lensed around a Schwarzschild black hole over the Gaia star field">
  <p>Warm disc (4500 K): the bright left rim is the Doppler-boosted approaching side; the thin arc over the shadow is the lensed far side plus the photon ring.</p>
</div>

<div align="center">
  <img src="schwarzschild-disc-stars-blue-1920x1080-2026-09-23.png" alt="A hot 15000 K blue-white blackbody accretion disc lensed around a Schwarzschild black hole over the Gaia star field">
  <p>Hot disc (15000 K): a blackbody this hot glows blue-white (it never reaches a saturated blue), so a brighter, more energetic-looking disc than the warm one.</p>
</div>

### Volumetric accretion disc over the Gaia star field

A thick (volumetric) blackbody accretion disc over the Gaia DR3 (G &le; 12) star
field, single pass. Unlike the razor-thin flat disc, the gas has real thickness,
so its edge is a soft, fluffy falloff (no sharp outer rim) and its opacity comes
from the local gas density: the dim receding side still fully occludes the stars
behind it (dim in emission, unchanged in opacity), while only the tenuous outer
gas lets starlight bleed through. The disc's `flux_scale` dims its emitted light
into the star field's range so one exposure holds both. Graded from linear HDR
with bloom + ACES (recipes:
`images/schwarzschild-voldisc-stars-yellow-1920x1080-2026-09-23.toml`,
`images/schwarzschild-voldisc-stars-blue-1920x1080-2026-09-23.toml`). An HDR
gain-map JPEG (`*-hdr.jpg`) sits beside each PNG: an ordinary SDR JPEG carrying a
hidden gain map, so it shows SDR everywhere and brightens the disc past white on
an HDR screen.

<div align="center">
  <img src="schwarzschild-voldisc-stars-yellow-1920x1080-2026-09-23.png" alt="A warm volumetric blackbody accretion disc around a Schwarzschild black hole over the Gaia star field">
  <p>Warm volumetric disc (4500 K). HDR version: <code>schwarzschild-voldisc-stars-yellow-1920x1080-2026-09-23-hdr.jpg</code>.</p>
</div>

<div align="center">
  <img src="schwarzschild-voldisc-stars-blue-1920x1080-2026-09-23.png" alt="A hot blue-white volumetric blackbody accretion disc around a Schwarzschild black hole over the Gaia star field">
  <p>Hot volumetric disc (15000 K, blue-white). HDR version: <code>schwarzschild-voldisc-stars-blue-1920x1080-2026-09-23-hdr.jpg</code>.</p>
</div>

### Checkerboard Accretion Disk
<div align="center">
  <img src="render_schwarzschild_checker_texture.png" alt="Schwarzschild black hole with a checkerboard accretion disk showing gravitational lensing effects">
  <p>Schwarzschild black hole with a checkerboard accretion disk showing gravitational lensing effects</p>
</div>

### Volumetric Rendering
<div align="center">
  <img src="render_schwarzschild_volumetric.png" alt="Schwarzschild black hole visualization with volumetric light disc effects">
  <p>Schwarzschild black hole visualization with volumetric disc (background image: <a href="https://commons.wikimedia.org/wiki/File:NGC6355_-_HST_-_Potw2301a.jpg">NGC 6355</a>)</p>
</div>

### Close vantages

First-person views from a static observer close to the hole (recipes:
`images/create-vantage-images.sh`; checker sky, 3500 K blackbody disc so
gravitational blueshift keeps the gas in warm hues). The [Kerr
gallery](kerr.md#close-vantages) has the same four placements around a
spinning hole for comparison.

<div align="center">
  <img src="vantage_porthole.png" alt="The whole sky compressed into a porthole, seen from just outside the horizon">
  <p>Porthole (r = 1.05 r_s, looking straight out): this close to the horizon the entire sky, disc included, is compressed into a shrinking circle of light; everything around it is the black hole filling the view.</p>
</div>

<div align="center">
  <img src="vantage_photonsphere_edge.png" alt="The winding wall of light seen from inside the photon sphere">
  <p>Photon-sphere edge (r = 1.25 r_s, looking sideways): near the photon sphere at 1.5 r_s, sideways-directed light circles the hole; the bright wall is made of ever-tighter wound images of the disc and sky.</p>
</div>

<div align="center">
  <img src="vantage_grazing_winding.png" alt="Multiple wound images of the sky in a grazing band">
  <p>Grazing view from just above the photon sphere (r = 1.55 r_s): the band stacks multiple complete windings of the sky, each successive image thinner than the last.</p>
</div>

<div align="center">
  <img src="vantage_rearview.png" alt="Looking away from the black hole: a magnified sky and the disc band">
  <p>Rearview (r = 2 r_s, looking away and slightly up): the outward sky, gravitationally magnified, with the accretion disc cutting through as a bright band.</p>
</div>
