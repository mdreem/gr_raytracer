//! Star tube subdivision, sky footprints, and flux gathering.

use super::progress::Progress;
use super::Raytracer;
use super::tube::SampleTube;
use crate::geometry::geometry::Geometry;
use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::CIETristimulus;
use crate::rendering::light::radiance::Radiance;
use crate::rendering::scene::{EscapeInfo, RayClass, RaySample};
use crate::rendering::stars::star_catalog::Star;
use log::{debug, warn};
use nalgebra::Vector3;
use rayon::prelude::*;
use std::ops::{Add, Mul};
use std::sync::atomic::{AtomicUsize, Ordering};

struct StarCollectionData {
    pub total_g_mag: f64,
}

impl Add for StarCollectionData {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        StarCollectionData {
            total_g_mag: self.total_g_mag + other.total_g_mag,
        }
    }
}

impl Mul<StarCollectionData> for f64 {
    type Output = StarCollectionData;

    fn mul(self, rhs: StarCollectionData) -> Self::Output {
        StarCollectionData {
            total_g_mag: self * rhs.total_g_mag,
        }
    }
}

fn solid_angle_from_vecs(
    a_vec: &Vector3<f64>,
    b_vec: &Vector3<f64>,
    c_vec: &Vector3<f64>,
) -> Result<f64, RaytracerError> {
    let num = a_vec.dot(&b_vec.cross(c_vec));
    let denom = 1.0 + a_vec.dot(b_vec) + b_vec.dot(c_vec) + c_vec.dot(a_vec);

    let angle = 2.0 * num.atan2(denom);
    // Check each triangle, not just the total quad area: one collapsed
    // triangle must not reach the star-membership test even if the other
    // has a valid area. Equal vertices can leave a tiny nonzero determinant
    // through cancellation, so check them explicitly too.
    if a_vec == b_vec
        || b_vec == c_vec
        || c_vec == a_vec
        || !num.is_finite()
        || num == 0.0
        || !denom.is_finite()
        || !angle.is_finite()
        || angle == 0.0
    {
        Err(RaytracerError::InvalidStarTube(
            "nonfinite or degenerate spherical triangle",
        ))
    } else {
        Ok(angle)
    }
}

pub(super) fn compute_traced_tube_solid_angle(
    a: &EscapeInfo,
    b: &EscapeInfo,
    c: &EscapeInfo,
    d: &EscapeInfo,
) -> Result<f64, RaytracerError> {
    let a_vec = Vector3::new(a.x, a.y, a.z);
    let b_vec = Vector3::new(b.x, b.y, b.z);
    let c_vec = Vector3::new(c.x, c.y, c.z);
    let d_vec = Vector3::new(d.x, d.y, d.z);

    compute_solid_angle(&a_vec, &b_vec, &c_vec, &d_vec)
}

fn compute_solid_angle(
    a_vec: &Vector3<f64>,
    b_vec: &Vector3<f64>,
    c_vec: &Vector3<f64>,
    d_vec: &Vector3<f64>,
) -> Result<f64, RaytracerError> {
    // See https://en.wikipedia.org/wiki/Solid_angle
    // |a_vec| = |b_vec| = |c_vec| = 1, so the formula simplifies to:

    // Triangulate the quad as tri(a,b,c) + tri(b,d,c): both traversed with the
    // same orientation (corners are TL,TR,BL,BR). Using (b,c,d) instead flips
    // the second triangle's orientation, so the two nearly cancel and the sum
    // is a twist residual rather than the quad area. Matches the gather's split.
    let angle_1 = solid_angle_from_vecs(a_vec, b_vec, c_vec)?;
    let angle_2 = solid_angle_from_vecs(b_vec, d_vec, c_vec)?;

    // Sum the absolute triangle areas: robust to a per-triangle sign flip at a
    // fold (where the signed sum would cancel), and gives the true magnitude.
    Ok(angle_1.abs() + angle_2.abs())
}

pub(super) fn finite_tube_flux(color: Vector3<f64>) -> Result<Vector3<f64>, RaytracerError> {
    if color.iter().all(|v| v.is_finite()) {
        Ok(color)
    } else {
        Err(RaytracerError::InvalidStarTube("nonfinite gathered color"))
    }
}

/// Angular flatness tolerance (radians) for refining a tube edge into a polyline
/// that follows the true sky arc: refine an edge segment while its traced
/// midpoint deviates from the straight chord by more than this.
const EDGE_FLATNESS_TOLERANCE: f64 = 0.003;

/// Recursion cap for edge refinement (a 1-D subdivision, so 2^depth points).
const MAX_EDGE_REFINE_DEPTH: u32 = 12;

impl<'a, G: Geometry> Raytracer<'a, G> {
    /// `depth` counts splits already made: zero is the base tube, and a
    /// configured maximum of one permits exactly one split.
    pub(super) fn subdivide(
        &self,
        sample_tube: &SampleTube<G::Chart>,
        depth: usize,
    ) -> Result<Option<Vector3<f64>>, RaytracerError> {
        if depth >= self.scene.max_subdivision_depth {
            debug!("Maximum subdivision depth reached.");
            return Ok(None);
        }

        let samples = self.trace_subdivision(sample_tube)?;
        let mut color = Vector3::zeros();
        for child in sample_tube.children(&samples) {
            color += self.compute_color(&child, depth + 1)?;
        }
        Ok(Some(finite_tube_flux(color)?))
    }

    pub(super) fn trace_subdivision(
        &self,
        sample_tube: &SampleTube<G::Chart>,
    ) -> Result<[RaySample<G::Chart>; 5], RaytracerError> {
        let points = sample_tube
            .subdivision_points()
            .ok_or(RaytracerError::InvalidStarTube(
                "tube screen bounds cannot be subdivided",
            ))?;
        let trace = |(row, col): (f64, f64)| {
            // The integer labels stay attached to the owning pixel; only
            // the tube's screen bounds determine the actual sampling position.
            let pixel_row = sample_tube.a.ray.row;
            let pixel_col = sample_tube.a.ray.col;
            let ray = self.scene.camera.get_ray_for_offset(
                pixel_row,
                pixel_col,
                col - pixel_col as f64 + 0.5,
                row - pixel_row as f64 + 0.5,
            );
            self.scene.color_of_ray(&ray)
        };
        // Five traces, with their results borrowed by all four children.
        Ok([
            trace(points[0])?,
            trace(points[1])?,
            trace(points[2])?,
            trace(points[3])?,
            trace(points[4])?,
        ])
    }

    pub(super) fn compute_color(
        &self,
        sample_tube: &SampleTube<G::Chart>,
        depth: usize,
    ) -> Result<Vector3<f64>, RaytracerError> {
        let color = match (
            sample_tube.a.ray_class,
            sample_tube.b.ray_class,
            sample_tube.c.ray_class,
            sample_tube.d.ray_class,
        ) {
            (
                RayClass::Escaped(a),
                RayClass::Escaped(b),
                RayClass::Escaped(c),
                RayClass::Escaped(d),
            ) => {
                let angles = [
                    sample_tube.a.accumulated_angular_distance,
                    sample_tube.b.accumulated_angular_distance,
                    sample_tube.c.accumulated_angular_distance,
                    sample_tube.d.accumulated_angular_distance,
                ];
                if !angles.iter().all(|angle| angle.is_finite()) {
                    return Err(RaytracerError::InvalidStarTube("nonfinite winding angle"));
                }
                let min_angle = angles.iter().copied().fold(f64::INFINITY, f64::min);
                let max_angle = angles.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let spread = max_angle - min_angle;
                if spread > self.scene.winding_spread_threshold
                    && let Some(color) = self.subdivide(sample_tube, depth)?
                {
                    return Ok(color);
                }
                let o_a = sample_tube
                    .a
                    .ray
                    .momentum
                    .to_cartesian()
                    .get_spatial_vector()
                    .normalize();
                let o_b = sample_tube
                    .b
                    .ray
                    .momentum
                    .to_cartesian()
                    .get_spatial_vector()
                    .normalize();
                let o_c = sample_tube
                    .c
                    .ray
                    .momentum
                    .to_cartesian()
                    .get_spatial_vector()
                    .normalize();
                let o_d = sample_tube
                    .d
                    .ray
                    .momentum
                    .to_cartesian()
                    .get_spatial_vector()
                    .normalize();

                if self.scene.curved_star_membership {
                    return finite_tube_flux(self.gather_curved(
                        sample_tube,
                        &a,
                        &b,
                        &c,
                        &d,
                        &o_a,
                        &o_b,
                        &o_c,
                        &o_d,
                        depth,
                    )?);
                }

                let original_angle = compute_solid_angle(&o_a, &o_b, &o_c, &o_d)?;
                let solid_angle = compute_traced_tube_solid_angle(&a, &b, &c, &d)?;
                let ratio = self.magnification(original_angle, solid_angle);
                if !ratio.is_finite() || ratio <= 0.0 {
                    return Err(RaytracerError::InvalidStarTube("invalid magnification"));
                }

                ratio * self.compute_star_collection_data(&a, &b, &c, &d)
            }
            // The star layer only carries background starlight. A tube with no
            // escaped corner sees no background (the hole or an opaque object
            // fills it), so it contributes zero star flux; the object/shadow
            // itself already comes from the base pass. Returning the corner's
            // foreground colour here instead would leak that colour (e.g. the
            // disc) into the star layer and paint a spurious rim along a
            // silhouette, resolution-dependently.
            (RayClass::Captured, RayClass::Captured, RayClass::Captured, RayClass::Captured)
            | (RayClass::Hit, RayClass::Hit, RayClass::Hit, RayClass::Hit) => Vector3::zeros(),
            // Mixed corners: subdivide so the escaped sub-tubes still gather
            // their stars; if subdivision is capped, contribute no star flux
            // rather than leaking foreground colour.
            _ => {
                let any_escaped = matches!(sample_tube.a.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.b.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.c.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.d.ray_class, RayClass::Escaped(_));
                if any_escaped {
                    self.subdivide(sample_tube, depth)?
                        .unwrap_or(Vector3::zeros())
                } else {
                    Vector3::zeros()
                }
            }
        };
        finite_tube_flux(color)
    }

    /// Lensing magnification `A/B` of a tube (screen solid angle over footprint
    /// solid angle), clamped to `star_magnification_cap`. The clamp bounds the
    /// caustic blow-up where `B` collapses toward zero.
    fn magnification(&self, original_angle: f64, footprint: f64) -> f64 {
        (original_angle / footprint).min(self.scene.star_magnification_cap)
    }

    fn compute_star_collection_data(
        &self,
        a: &EscapeInfo,
        b: &EscapeInfo,
        c: &EscapeInfo,
        d: &EscapeInfo,
    ) -> Vector3<f64> {
        // Accumulate each in-tube star's XYZ radiance (flux times its blackbody
        // chromaticity), so hues sum in linear light. Y carries the summed
        // flux; X and Z carry the colour.
        // Flux has no opacity: magnification and subdivision must scale/sum
        // light only, never an alpha that later weights that same light again.
        let mut total = Vector3::zeros();
        // One frequency shift per tube: the four corners share nearly the same
        // sky direction, so average their g. (Away from the photon ring g is
        // near-uniform; the winding subdivision already splits the tubes where
        // it is not.)
        let redshift = 0.25 * (a.redshift + b.redshift + c.redshift + d.redshift);
        // TODO: Move the collection into scene
        if let Some(star_catalog) = &self.scene.star_catalog {
            // The two triangles tile the tube's quad. The octree returns exactly
            // the stars inside each (its plane test is the spherical-triangle
            // test), so a star shared on the diagonal is gathered by both, which
            // reproduces the old double-hit weighting.
            let triangles = [
                [a.to_vec(), b.to_vec(), c.to_vec()],
                [b.to_vec(), d.to_vec(), c.to_vec()],
            ];
            let mut stars = Vec::new();
            for triangle in &triangles {
                stars.clear();
                star_catalog.stars.gather_triangle(triangle, &mut stars);
                for star in &stars {
                    let emission = self.redshifted_star_emission(star, redshift);
                    total.x += emission.x;
                    total.y += emission.y;
                    total.z += emission.z;
                }
            }
        }

        let scale = self.scene.star_flux_scale;
        total * scale
    }

    /// Curved-boundary star gather for an all-escaped tube. Splits a folded tube
    /// (the critical curve passes through it) so each piece is single-sheet,
    /// then gathers against a boundary polygon whose edges follow the true sky
    /// arcs instead of straight chords. One gather per tube (no re-gather into
    /// sub-tubes), so it fills ring gaps without concentrating or thinning flux.
    #[allow(clippy::too_many_arguments)]
    fn gather_curved(
        &self,
        tube: &SampleTube<G::Chart>,
        a: &EscapeInfo,
        b: &EscapeInfo,
        c: &EscapeInfo,
        d: &EscapeInfo,
        o_a: &Vector3<f64>,
        o_b: &Vector3<f64>,
        o_c: &Vector3<f64>,
        o_d: &Vector3<f64>,
        depth: usize,
    ) -> Result<Vector3<f64>, RaytracerError> {
        let (va, vb, vc, vd) = (a.to_vec(), b.to_vec(), c.to_vec(), d.to_vec());
        // Fold test: the two tiling triangles wind opposite ways when the
        // critical curve crosses this tube. Split (reusing the subdivider, whose
        // children re-enter this path) until each piece is unfolded or capped.
        let folded = va.dot(&vb.cross(&vc)) * vb.dot(&vd.cross(&vc)) < 0.0;
        let original_angle = compute_solid_angle(o_a, o_b, o_c, o_d)?;
        if folded {
            // The critical curve crosses this tube. Split it so each piece is a
            // single sheet; but subdivision can still drop a star into a sliver
            // between child polygons, so never gather less than the flat path:
            // take whichever of the split result and the flat gather is larger.
            if depth < self.scene.max_subdivision_depth
                && let Some(color) = self.subdivide(tube, depth)?
            {
                let flat = self.flat_gather(a, b, c, d, original_angle)?;
                return Ok(if color.norm() >= flat.norm() {
                    color
                } else {
                    flat
                });
            }
            // Could not split (depth cap): the curved polygon would
            // self-intersect and its fan gather can drop stars, so use flat.
            return self.flat_gather(a, b, c, d, original_angle);
        }
        match self.tube_boundary(tube, &va, &vb, &vd, &vc) {
            Some(boundary) if boundary.len() >= 3 => {
                let redshift = 0.25 * (a.redshift + b.redshift + c.redshift + d.redshift);
                let (flux, footprint, signed_area) = self.gather_polygon(&boundary, redshift);
                // Hidden fold: an even number of critical-curve crossings the
                // corner parity test misses, so the polygon still self-intersects
                // and the signed area cancels. Fall back to the flat gather.
                if footprint > 0.0 && signed_area.abs() < 0.5 * footprint {
                    return self.flat_gather(a, b, c, d, original_angle);
                }
                if footprint > 0.0 && footprint.is_finite() {
                    return Ok(self.magnification(original_angle, footprint) * flux);
                }
                Ok(Vector3::zeros())
            }
            // An edge crossed the shadow (a corner failed to escape) or the
            // bounds could not be refined: fall back to the flat-quad gather.
            _ => self.flat_gather(a, b, c, d, original_angle),
        }
    }

    /// Trace a single screen point of the owning pixel; Some(dir) if it escapes.
    fn trace_escape(
        &self,
        pixel_row: i64,
        pixel_col: i64,
        row: f64,
        col: f64,
    ) -> Option<Vector3<f64>> {
        let ray = self.scene.camera.get_ray_for_offset(
            pixel_row,
            pixel_col,
            col - pixel_col as f64 + 0.5,
            row - pixel_row as f64 + 0.5,
        );
        match self.scene.color_of_ray(&ray).ok()?.ray_class {
            RayClass::Escaped(e) => Some(e.to_vec()),
            _ => None,
        }
    }

    /// Points strictly after `d0` up to and including `d1`, following the true
    /// arc: recurse while the traced midpoint deviates from the chord.
    fn refine_edge(
        &self,
        pr: i64,
        pc: i64,
        p0: (f64, f64),
        d0: Vector3<f64>,
        p1: (f64, f64),
        d1: Vector3<f64>,
        depth: u32,
    ) -> Option<Vec<Vector3<f64>>> {
        if depth >= MAX_EDGE_REFINE_DEPTH {
            return Some(vec![d1]);
        }
        let pm = ((p0.0 + p1.0) * 0.5, (p0.1 + p1.1) * 0.5);
        let dm = self.trace_escape(pr, pc, pm.0, pm.1)?;
        let chord_mid = (d0 + d1).normalize();
        let deviation = dm.dot(&chord_mid).clamp(-1.0, 1.0).acos();
        if deviation <= EDGE_FLATNESS_TOLERANCE {
            return Some(vec![d1]);
        }
        let mut left = self.refine_edge(pr, pc, p0, d0, pm, dm, depth + 1)?;
        left.extend(self.refine_edge(pr, pc, pm, dm, p1, d1, depth + 1)?);
        Some(left)
    }

    /// Refined boundary polygon of a tube (corners in order TL, TR, BR, BL, i.e.
    /// `va, vb, vd, vc`), each edge a polyline hugging the true sky arc.
    fn tube_boundary(
        &self,
        tube: &SampleTube<G::Chart>,
        va: &Vector3<f64>,
        vb: &Vector3<f64>,
        vd: &Vector3<f64>,
        vc: &Vector3<f64>,
    ) -> Option<Vec<Vector3<f64>>> {
        let sb = tube.screen_bounds;
        let pr = tube.a.ray.row;
        let pc = tube.a.ray.col;
        let corners = [
            ((sb.top, sb.left), *va),
            ((sb.top, sb.right), *vb),
            ((sb.bottom, sb.right), *vd),
            ((sb.bottom, sb.left), *vc),
        ];
        let mut poly = vec![corners[0].1];
        for i in 0..4 {
            let (p0, d0) = corners[i];
            let (p1, d1) = corners[(i + 1) % 4];
            poly.extend(self.refine_edge(pr, pc, p0, d0, p1, d1, 0)?);
        }
        poly.pop(); // drop the closing duplicate of corner a
        Some(poly)
    }

    /// Fan-triangulate the boundary polygon from its centroid; sum star flux
    /// (each star once) and the enclosed solid angle (the true footprint area).
    /// Returns `(flux, footprint, signed_area)`: the summed star flux, the
    /// footprint solid angle (sum of |triangle areas|), and the signed area sum.
    /// A simple polygon has `|signed_area| == footprint`; a self-intersecting
    /// (folded) one has the signs cancel, so `|signed_area| << footprint` flags
    /// a fold the fan gather cannot be trusted on.
    fn gather_polygon(&self, poly: &[Vector3<f64>], redshift: f64) -> (Vector3<f64>, f64, f64) {
        let n = poly.len();
        let mut centroid = Vector3::zeros();
        for p in poly {
            centroid += p;
        }
        let centroid = centroid.normalize();
        let mut flux = Vector3::zeros();
        let mut footprint = 0.0;
        let mut signed_area = 0.0;
        if let Some(catalog) = &self.scene.star_catalog {
            let mut seen = std::collections::HashSet::new();
            let mut stars = Vec::new();
            for i in 0..n {
                let p0 = poly[i];
                let p1 = poly[(i + 1) % n];
                if let Ok(area) = solid_angle_from_vecs(&centroid, &p0, &p1) {
                    footprint += area.abs();
                    signed_area += area;
                }
                stars.clear();
                catalog
                    .stars
                    .gather_triangle(&[centroid, p0, p1], &mut stars);
                for star in &stars {
                    if seen.insert(star.source_id) {
                        let emission = self.redshifted_star_emission(star, redshift);
                        flux.x += emission.x;
                        flux.y += emission.y;
                        flux.z += emission.z;
                    }
                }
            }
        }
        (flux * self.scene.star_flux_scale, footprint, signed_area)
    }

    /// The original flat-quad gather: the star flux inside the tube's four
    /// corners, magnified by the corner solid-angle ratio. Used as the fallback
    /// for tubes whose fold cannot be resolved, so the curved path never gathers
    /// less than the flat path.
    fn flat_gather(
        &self,
        a: &EscapeInfo,
        b: &EscapeInfo,
        c: &EscapeInfo,
        d: &EscapeInfo,
        original_angle: f64,
    ) -> Result<Vector3<f64>, RaytracerError> {
        let solid_angle = compute_traced_tube_solid_angle(a, b, c, d)?;
        let ratio = self.magnification(original_angle, solid_angle);
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(RaytracerError::InvalidStarTube("invalid magnification"));
        }
        Ok(ratio * self.compute_star_collection_data(a, b, c, d))
    }

    /// One star's observed XYZ radiance under the frequency shift `g`.
    ///
    /// A blackbody at `T` seen with shift `g` is a blackbody at `g * T` (Wien),
    /// so the hue comes from the LUT at `g * T`, Y-normalised. Brightness comes
    /// from the catalogue flux boosted by `g^4` (the bolometric surface-
    /// brightness law `I_obs = g^4 I_emit`); the solid-angle magnification is
    /// applied separately by the tube `ratio`, so it is not double-counted
    /// here. `g == 1` (or no LUT) returns the precomputed rest-frame emission.
    fn redshifted_star_emission(&self, star: &Star, g: f64) -> CIETristimulus {
        if (g - 1.0).abs() < 1e-6 {
            return star.emission_xyz;
        }
        let Some(mapper) = &self.scene.star_blackbody else {
            return star.emission_xyz;
        };
        let bb = mapper.blackbody_xyz(g * star.temperature, 1.0);
        let inv_y = if bb.y > 0.0 { 1.0 / bb.y } else { 0.0 };
        let luminance = star.relative_flux() * g.powi(4);
        CIETristimulus::new(
            bb.x * inv_y * luminance,
            luminance,
            bb.z * inv_y * luminance,
            1.0,
        )
    }

    pub(super) fn handle_tube(
        &self,
        sample_tube: &SampleTube<G::Chart>,
    ) -> Result<Vector3<f64>, RaytracerError> {
        self.compute_color(sample_tube, 0)
    }

    pub(super) fn add_star_layer(
        &self,
        buffer: &[RaySample<G::Chart>],
        width: usize,
        _height: usize,
        colors: &mut [Radiance],
    ) -> Result<(), RaytracerError> {
        if self
            .scene
            .star_catalog
            .as_ref()
            .is_none_or(|catalog| catalog.is_empty())
        {
            return Ok(());
        }
        // Retain the current center-to-center base footprints. Moving them
        // to pixel boundaries is separate from fixing recursive subdivision.
        // Each tube writes only its own pixel and reads the shared buffer, so
        // the gather parallelizes over pixels like the geodesic pass. This
        // matters most for the curved-membership path, whose per-tube cost near
        // the critical curve would otherwise pin one core.
        let tubes = AtomicUsize::new(0);
        let skipped = AtomicUsize::new(0);
        let gathered = AtomicUsize::new(0);
        let progress = Progress::new(colors.len() as u64, "stars")?;
        colors.par_iter_mut().enumerate().for_each(|(idx, color)| {
            progress.set(gathered.fetch_add(1, Ordering::Relaxed) as u64 + 1);
            let (row, col) = (idx / width, idx % width);
            // The base tube spans (row, col)..(row+1, col+1); on the last row or
            // column `from_buffer` clamps the stencil inward so the boundary
            // pixel still gets its inward neighbour's footprint instead of being
            // skipped (which would drop border stars and leave section seams).
            if let Some(tube) =
                SampleTube::from_buffer(buffer, row as u32, col as u32, width as u32)
            {
                tubes.fetch_add(1, Ordering::Relaxed);
                // Skip a tube whose gather fails rather than aborting the whole
                // section; it just contributes no star flux.
                match self.handle_tube(&tube) {
                    Ok(stars) => {
                        // Star layer goes UNDER the foreground, weighted by the
                        // foreground's surviving transmittance.
                        *color = color.over(Radiance::new(stars.x, stars.y, stars.z, 0.0));
                    }
                    Err(err) => {
                        skipped.fetch_add(1, Ordering::Relaxed);
                        debug!("Skipping star gather for {:?}: {}", tube.screen_bounds, err);
                    }
                }
            }
        });
        progress.finish();
        let tubes = tubes.into_inner();
        let skipped = skipped.into_inner();
        if skipped > 0 {
            warn!(
                "Star gather skipped {} of {} tubes ({:.2}%) in this section; \
                 those tubes contributed no star flux",
                skipped,
                tubes,
                100.0 * skipped as f64 / tubes as f64
            );
        }
        Ok(())
    }
}
