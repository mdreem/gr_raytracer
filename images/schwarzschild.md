# Schwarzschild Black Hole

[← Back to the gallery index](images.md)

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

### Gaia star field

A Schwarzschild black hole lensing the real Gaia DR3 catalogue (G &le; 12) as
point-source stars, no accretion disc. The shadow and photon ring are perfectly
circular; the ring glow is clumpy only because the real sky is uneven (recipe:
`images/gaia-mag12-schwarzschild-no-disc-1920x1080-2026-09-21.toml`; rendered to
linear HDR and graded with bloom + an ACES tone map at exposure 69).

<div align="center">
  <img src="gaia-mag12-schwarzschild-no-disc-1920x1080-2026-09-21.png" alt="A Schwarzschild black hole gravitationally lensing the real Gaia DR3 star field">
  <p>Gaia DR3 stars (G &le; 12) lensed by a Schwarzschild black hole; 1920&times;1080 HDR, bloom + ACES.</p>
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
