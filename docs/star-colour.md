# Star colour (blackbody tinting of catalogue point sources)

How a catalogue star's Gaia photometry becomes an on-screen colour, why it is
done this way, and the more accurate method we can swap in later.

Relevant code: `src/rendering/star_catalog.rs` (per-star precompute) and
`src/rendering/raytracer.rs::compute_star_collection_data` (the gather).

## The governing principle: temperature sets colour, not brightness

The accretion disc uses `black_body_radiation::integrate_blackbody_xyz(T, z)`
at its **absolute** amplitude: a hotter annulus is genuinely brighter, because
the disc's brightness *is* its surface temperature (Planck radiance climbs
steeply with `T`).

A star is different. Its on-screen brightness is set by its **observed
magnitude** (luminosity / distance^2), not by its surface temperature. A
distant hot O star and a nearby cool M dwarf can share the same observed flux.
So a star must take:

- its **brightness** from the catalogue magnitude (`flux = 10^(-0.4 * G)`), and
- only its **colour (chromaticity)** from the temperature.

Concretely (`star_catalog::star_emission_xyz`):

1. `bp_rp -> T` (colour temperature, see below).
2. `bb = get_cie_xyz_of_black_body_redshifted(T, 1.0)` -- absolute Planck XYZ.
3. Divide by `bb.y` to strip the absolute amplitude and keep only the hue
   (a Y-normalised chromaticity: `Y = 1`, `X`/`Z` carry the colour).
4. Multiply by the Pogson flux, giving a per-star XYZ whose luminance `Y`
   equals the observed flux.

The gather sums these XYZ vectors (linear light adds), then applies the global
`flux_scale` and the achromatic lensing magnification `ratio`. Magnification
scales the vector without changing its hue, which is correct: gravitational
lensing is achromatic.

Because the emission depends only on `g_mag` and `bp_rp`, it is precomputed
once per star at load time (`Star::emission_xyz`).

## `bp_rp -> T`: the polynomial (current implementation)

`colour_temperature_from_bp_rp` uses a Ballesteros-style reciprocal:

```
T = a / (bp_rp + b) + c,   a = 8235, b = 0.997, c = 1240   (Kelvin)
```

Why a reciprocal and not a literal polynomial in `bp_rp`: by Wien's law the
colour index runs roughly as `1/T`, so a reciprocal is the natural form and
stays well-behaved at the blue and red ends, where a polynomial fit would
over/undershoot (even go negative).

The three constants are pinned to reference `(BP-RP, Teff)` anchors spanning
the sequence, from the Pecaut & Mamajek dwarf colour-temperature table plus the
solar point:

| BP-RP | Teff (K) | note        |
|-------|----------|-------------|
| 0.00  | 9500     | hot end (A) |
| 0.82  | 5772     | Sun (G2V)   |
| 3.00  | 3300     | cool end (M)|

It reproduces those exactly and stays within a few percent across O..M, which
is ample for a display colour. It is **not** a precision Teff estimator (no
metallicity, no giant/dwarf distinction, no reddening correction). The input is
clamped to `bp_rp in [-0.5, 5.0]` so a stray very-blue value cannot drive the
denominator toward zero.

## `bp_rp -> T`: colour temperature (more accurate, deferred)

The polynomial gives the star's approximate *physical* effective temperature.
But we render each star **as a blackbody**, and a real star's `bp_rp` differs
from a blackbody's at the same physical `Teff` (line blanketing, the Balmer
jump). Feeding the physical `Teff` into a blackbody therefore yields a subtly
wrong displayed colour.

The self-consistent fix is the **colour temperature**: the temperature of the
blackbody whose *own synthetic Gaia `bp_rp`* matches the observed `bp_rp`. A
blackbody rendered at that temperature reproduces the observed colour index by
construction.

What it needs beyond the polynomial:

1. **The Gaia BP and RP passband curves** -- transmission vs wavelength,
   `[(lambda, S_BP, S_RP)]`, ~330..1050 nm (RP runs redder than the 830 nm the
   CMF integrator stops at, so this is a separate, wider integration). Source:
   SVO Filter Profile Service (`GAIA/GAIA3.Gbp`, `GAIA/GAIA3.Grp`) or the
   official "Gaia DR3 passbands" table. A few hundred rows; embed or drop in
   `resources/`.
2. **One zero-point constant** `C`. Because `bp_rp = BP_mag - RP_mag` and each
   magnitude has its own zero-point, matching the observed scale needs the
   BP-RP zero-point offset. Absorb it by anchoring: choose `C` so the Sun
   (`T = 5772 K`) reproduces `bp_rp = 0.82`. That single calibration also soaks
   up the photon-vs-energy detail below.

Sketch:

```rust
// Gaia is a photon counter: weight the energy radiance by lambda (∝ photon
// rate). hc cancels in the ratio; a uniform-grid dλ cancels too.
fn synth_bp_rp(t: f64) -> f64 {
    let (mut n_bp, mut n_rp) = (0.0, 0.0);
    for &(lambda_nm, s_bp, s_rp) in PASSBANDS {          // the extra data
        let lambda = lambda_nm * 1e-9;
        let photons = planck_spectral_radiance(lambda, t) * lambda;
        n_bp += photons * s_bp;
        n_rp += photons * s_rp;
    }
    -2.5 * (n_bp / n_rp).log10() + C                     // C from solar anchor
}

// Build once: monotone (cooler -> redder -> larger bp_rp), then invert by
// binary search + lerp, clamped to the grid ends.
fn build_lut() -> Vec<(f64 /*bp_rp*/, f64 /*T*/)> { /* over a log-T grid */ }
```

`planck_spectral_radiance` is currently private in `black_body_radiation.rs`;
exposing it (`pub(crate)`) plus adding the passband table is the only plumbing.
`colour_temperature_from_bp_rp` is the single swap point -- nothing at the call
sites changes.

## Redshift interaction (separate future item)

The black-hole redshift makes `T_obs = z * T`, which shifts the hue. When that
lands, look up the chroma at `z * T` per tube (the star keeps its rest-frame
`temperature`, which is why `Star::temperature` is retained) instead of using
the precomputed `emission_xyz`. The `z^5` amplitude boost that
`integrate_blackbody_xyz` bakes in should be stripped (we normalise by `Y`
anyway) and the flux boost applied separately as the `g`-factor on the flux, so
brightness stays governed by the catalogue, not the Planck amplitude.

## Rigour ceiling

Setting `Y = flux` treats the Gaia G-band flux as photopic luminance. They
differ by a temperature-dependent offset (G band != the eye's V response) of
order a few tenths of a magnitude. The exact version anchors the blackbody's
absolute scale so its *synthetic G-band* flux equals `10^(-0.4 G)`, then
integrates to XYZ for a true photopic `Y`. Invisible for a visualisation, but
noted as the fully rigorous form.
