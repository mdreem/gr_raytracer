use crate::configuration::AdaptiveSamplingConfig;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::CameraError;
use crate::rendering::color::{
    CIETristimulus, ToneMappingMethod, linear_srgb_to_srgb_buffer, xyz_to_linear_srgb,
    xyz_to_linear_srgb_buffer,
};
use crate::rendering::integrator::{IntegrationError, StopReason};
use crate::rendering::radiance::{Radiance, RadianceMean};
use crate::rendering::ray::{IntegratedRay, Ray};
use crate::rendering::scene::{EscapeInfo, RayClass, RaySample, Scene};
use crate::rendering::star_catalog::Star;
use crate::rendering::texture::TextureError;
use crate::rendering::tubetracer::SampleTube;
use image::{ImageBuffer, ImageError, ImageFormat, Rgb};
use indicatif::style::TemplateError;
use log::{debug, error, info, warn};
use nalgebra::Vector3;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::IntoParallelRefMutIterator;
use rayon::iter::ParallelIterator;
use std::io;
use std::ops::{Add, Mul};
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Debug, thiserror::Error)]
pub enum RaytracerError {
    #[error("Integration error: {0}")]
    IntegrationError(#[from] IntegrationError),
    #[error("I/O error: {0}")]
    IoError(#[from] io::Error),
    #[error("Configuration file error: {0}")]
    ConfigurationFileError(io::Error),
    #[error("TOML error: {0}")]
    TomlError(#[from] toml::de::Error),
    #[error("Texture error: {0}")]
    TextureError(#[from] TextureError),
    #[error("Unphysical nonpositive redshift factor {0} in beaming")]
    UnphysicalRedshift(f64),
    #[error("Image error: {0}")]
    ImageError(#[from] ImageError),
    #[error("Image buffer creation failed")]
    ImageBufferCreation,
    #[error("Progress bar template error: {0}")]
    ProgressBarTemplateError(#[from] TemplateError),
    #[error("Camera error: {0}")]
    CameraError(#[from] CameraError),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
    #[error("No circular orbit possible")]
    NoCircularOrbitPossible,
    #[error("Radius is below RISCO")]
    BelowRISCO,
    #[error("Radius is not finite")]
    NonFiniteRadius,
    #[error("Number is below zero")]
    NumberBelowZero,
    #[error("Denominator is close to zero")]
    DenominatorCloseToZero,
    #[error("Invalid star tube: {0}")]
    InvalidStarTube(&'static str),
}

pub struct Raytracer<'a, G: Geometry> {
    pub scene: Scene<'a, G>,
    tone_mapping: ToneMappingMethod,
    exposure: f64,
}

const MICHELSON_DENOMINATOR_EPSILON: f64 = 1e-4;

#[derive(PartialEq)]
struct PixelToSample {
    pub row: u32,
    pub col: u32,
    pub result: Option<Radiance>,
}

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

/// Michelson (relative) luminance contrast; relative because HDR radiance is unbounded.
fn luminance_contrast(p: &Radiance, q: &Radiance) -> f64 {
    let l_p = p.y;
    let l_q = q.y;
    (l_p - l_q).abs() / (l_p + l_q + MICHELSON_DENOMINATOR_EPSILON)
}

/// Absolute opacity difference equals the absolute transmittance difference.
fn opacity_contrast(p: &Radiance, q: &Radiance) -> f64 {
    (p.transmittance - q.transmittance).abs()
}

/// Faintness gate for the contrast triggers: true if the brighter of the pair
/// clears the floor (max keeps it symmetric). The floor is in linear CIE Y.
/// Scenes park the emissive disc far above unit luminance (Y > 1e4), while the
/// star-field background sits below ~1, so a floor of 1.0 keeps the disc fully
/// but stops the high-frequency background from flooding the contrast trigger
/// (it would otherwise flag ~all pixels). Class-change edges are not gated, so
/// silhouettes and the shadow rim are unaffected.
fn visible(p: &Radiance, q: &Radiance, minimum_luminance: f64) -> bool {
    p.y.max(q.y) > minimum_luminance
}

fn should_supersample_pair(
    pixel: &RaySample,
    neighbor: &RaySample,
    config: &AdaptiveSamplingConfig,
    minimum_luminance: f64,
) -> bool {
    if std::mem::discriminant(&pixel.ray_class) != std::mem::discriminant(&neighbor.ray_class) {
        return true;
    }

    if config.exclude_background_contrast && matches!(pixel.ray_class, RayClass::Escaped(_)) {
        return false;
    }

    visible(&pixel.color, &neighbor.color, minimum_luminance)
        && (luminance_contrast(&pixel.color, &neighbor.color) > config.luminance_contrast_threshold
            || opacity_contrast(&pixel.color, &neighbor.color) > config.opacity_contrast_threshold)
}

/// Fraction of the frame's 99th-percentile luminance used as the contrast
/// faintness floor when `minimum_luminance` is not set explicitly.
const RELATIVE_MINIMUM_LUMINANCE_FRACTION: f64 = 1e-3;

/// Resolve the contrast faintness floor: the configured absolute value, or a
/// scene-relative one. The relative floor is a small fraction of the
/// 99th-percentile luminance, which tracks the disc brightness while ignoring
/// firefly outliers (using `max` would be dominated by them).
fn resolve_minimum_luminance(config: &AdaptiveSamplingConfig, buffer: &[RaySample]) -> f64 {
    if let Some(value) = config.minimum_luminance {
        return value;
    }
    if buffer.is_empty() {
        return 0.0;
    }
    let mut luminances: Vec<f64> = buffer.iter().map(|sample| sample.color.y).collect();
    let index = (((luminances.len() - 1) as f64) * 0.99) as usize;
    luminances.select_nth_unstable_by(index, f64::total_cmp);
    RELATIVE_MINIMUM_LUMINANCE_FRACTION * luminances[index]
}

// https://rosettacode.org/wiki/Pseudo-random_numbers/Splitmix64
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn hash_pixel_samples(row: i64, col: i64, k: usize) -> f64 {
    let z = mix64((row as u64).wrapping_add(mix64((col as u64).wrapping_add(mix64(k as u64)))));

    // Convert to a floating-point number in the range [0, 1)
    (z >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

fn stratified_sample_offset(
    row: i64,
    col: i64,
    stratum_row: usize,
    stratum_col: usize,
    samples_per_axis: usize,
) -> (f64, f64) {
    debug_assert!(samples_per_axis > 0);
    let sample_index = stratum_row * samples_per_axis + stratum_col;
    let dx = (stratum_col as f64 + hash_pixel_samples(row, col, 2 * sample_index))
        / samples_per_axis as f64;
    let dy = (stratum_row as f64 + hash_pixel_samples(row, col, 2 * sample_index + 1))
        / samples_per_axis as f64;
    (dx, dy)
}

fn to_cartesian(theta: f64, phi: f64) -> nalgebra::Vector3<f64> {
    let x = theta.sin() * phi.cos();
    let y = theta.sin() * phi.sin();
    let z = theta.cos();
    nalgebra::Vector3::new(x, y, z)
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

fn compute_traced_tube_solid_angle(
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

fn finite_tube_flux(color: Vector3<f64>) -> Result<Vector3<f64>, RaytracerError> {
    if color.iter().all(|v| v.is_finite()) {
        Ok(color)
    } else {
        Err(RaytracerError::InvalidStarTube("nonfinite gathered color"))
    }
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    pub fn new(scene: Scene<'a, G>, tone_mapping: ToneMappingMethod, exposure: f64) -> Self {
        Self {
            scene,
            tone_mapping,
            exposure,
        }
    }

    #[allow(dead_code)] // For testing
    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        debug!("ray: {:?}", ray);
        let sample = self.scene.color_of_ray(&ray);
        debug!("sample: {:?}", sample);
    }

    /// `depth` counts splits already made: zero is the base tube, and a
    /// configured maximum of one permits exactly one split.
    fn subdivide(
        &self,
        sample_tube: &SampleTube,
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

    fn trace_subdivision(
        &self,
        sample_tube: &SampleTube,
    ) -> Result<[RaySample; 5], RaytracerError> {
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

    fn compute_color(
        &self,
        sample_tube: &SampleTube,
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
                    .get_cartesian_vector(&sample_tube.a.ray.position)
                    .normalize();
                let o_b = sample_tube
                    .b
                    .ray
                    .momentum
                    .get_cartesian_vector(&sample_tube.b.ray.position)
                    .normalize();
                let o_c = sample_tube
                    .c
                    .ray
                    .momentum
                    .get_cartesian_vector(&sample_tube.c.ray.position)
                    .normalize();
                let o_d = sample_tube
                    .d
                    .ray
                    .momentum
                    .get_cartesian_vector(&sample_tube.d.ray.position)
                    .normalize();

                let original_angle = compute_solid_angle(&o_a, &o_b, &o_c, &o_d)?;
                let solid_angle = compute_traced_tube_solid_angle(&a, &b, &c, &d)?;
                let ratio = original_angle / solid_angle;
                if !ratio.is_finite() || ratio <= 0.0 {
                    return Err(RaytracerError::InvalidStarTube("invalid magnification"));
                }

                ratio * self.compute_star_collection_data(&a, &b, &c, &d)
            }
            // TODO: These branches use corner foreground light as a star-flux
            // fallback. Correct composition does not make this classification
            // physically correct.
            // all four corners captured or all four hitting an opaque object:
            // Just return the color of corner a, no need to subdivide.
            (RayClass::Captured, RayClass::Captured, RayClass::Captured, RayClass::Captured)
            | (RayClass::Hit, RayClass::Hit, RayClass::Hit, RayClass::Hit) => {
                sample_tube.a.color.as_vector()
            }
            // Catch mixed cases.
            // If mixed but not-escaped: return the color of corner a, no need to subdivide.
            // If any corner escaped: subdivide. Of maximum is reached, just return the color of corner a.
            _ => {
                let any_escaped = matches!(sample_tube.a.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.b.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.c.ray_class, RayClass::Escaped(_))
                    || matches!(sample_tube.d.ray_class, RayClass::Escaped(_));
                if any_escaped {
                    self.subdivide(sample_tube, depth)?
                        .unwrap_or(sample_tube.a.color.as_vector())
                } else {
                    sample_tube.a.color.as_vector()
                }
            }
        };
        finite_tube_flux(color)
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

    fn handle_tube(&self, sample_tube: &SampleTube) -> Result<Vector3<f64>, RaytracerError> {
        self.compute_color(sample_tube, 0)
    }

    fn add_star_layer(
        &self,
        buffer: &[RaySample],
        width: usize,
        height: usize,
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
        let mut tubes = 0usize;
        let mut skipped = 0usize;
        for row in 0..height.saturating_sub(1) {
            for col in 0..width.saturating_sub(1) {
                if let Some(tube) =
                    SampleTube::from_buffer(buffer, row as u32, col as u32, width as u32)
                {
                    tubes += 1;
                    let idx = row * width + col;
                    // Skip a tube whose gather fails rather than aborting the whole
                    // section; it just contributes no star flux.
                    match self.handle_tube(&tube) {
                        Ok(stars) => {
                            // Star layer goes UNDER the foreground, weighted by
                            // the foreground's surviving transmittance.
                            colors[idx] =
                                colors[idx].over(Radiance::new(stars.x, stars.y, stars.z, 0.0));
                        }
                        Err(err) => {
                            skipped += 1;
                            debug!(
                                "Skipping star gather for {:?}: {}",
                                tube.screen_bounds, err
                            );
                        }
                    }
                }
            }
        }
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

    fn render_section_to_radiance_buffer(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<Radiance>, RaytracerError> {
        let buffer =
            self.render_section_to_radiance_buffer_raw(from_row, from_col, to_row, to_col)?;
        if self.scene.adaptive_sampling.enabled || self.scene.sampling_mask_color.is_some() {
            self.render_section_to_radiance_buffer_supersampled(
                from_row, from_col, to_row, to_col, &buffer,
            )
        } else {
            let width = (to_col - from_col) as usize;
            let height = (to_row - from_row) as usize;
            let mut colors: Vec<Radiance> = buffer.iter().map(|s| s.color).collect();
            self.add_star_layer(&buffer, width, height, &mut colors)?;

            Ok(colors)
        }
    }

    fn render_section_to_radiance_buffer_raw(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
    ) -> Result<Vec<RaySample>, RaytracerError> {
        let count = AtomicUsize::new(0);
        let max_count = (to_row - from_row) * (to_col - from_col);
        // Create some dummy rays. They will be replaced with the actual rays in the parallel loop below.
        let mut buffer: Vec<RaySample> = vec![
            RaySample {
                color: Radiance::BLACK,
                ray_class: RayClass::Escaped(EscapeInfo {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    redshift: 1.0,
                }),
                accumulated_angular_distance: 0.0,
                ray: Ray::new(
                    0,
                    0,
                    Point::new(
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        crate::geometry::point::CoordinateSystem::Cartesian,
                    ),
                    FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0),
                ),
            };
            max_count as usize
        ];

        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(max_count as u64);
        pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
            .map_err(RaytracerError::ProgressBarTemplateError)?
            .progress_chars("█▇▆▅▄▃▂▁  "));

        buffer.par_iter_mut().enumerate().for_each(|(i, p)| {
            count.fetch_add(1, Ordering::SeqCst);
            pb.set_position(count.load(Ordering::Relaxed) as u64);

            let y = i as u32 / (to_col - from_col);
            let x = i as u32 % (to_col - from_col);

            let ray = self
                .scene
                .camera
                .get_ray_for((y + from_row) as i64, (x + from_col) as i64);

            match self.scene.color_of_ray(&ray) {
                Ok(sample) => *p = sample,
                Err(err) => {
                    error!(
                        "Unable to compute color for ray at pixel ({}, {}): {:?}",
                        x + from_col,
                        y + from_row,
                        err
                    );
                }
            }
        });
        pb.finish();
        Ok(buffer)
    }

    fn get_pixel_index(
        &self,
        row: u32,
        col: u32,
        width: u32,
        offset_row: u32,
        offset_col: u32,
    ) -> usize {
        ((row - offset_row) * width + (col - offset_col)) as usize
    }

    fn render_section_to_radiance_buffer_supersampled(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer: &[RaySample],
    ) -> Result<Vec<Radiance>, RaytracerError> {
        info!(
            "Rendering section from ({}, {}) to ({}, {}) with supersampling",
            from_row, from_col, to_row, to_col
        );
        let samples_per_axis = self.scene.adaptive_sampling.samples_per_axis;
        let minimum_luminance = resolve_minimum_luminance(&self.scene.adaptive_sampling, buffer);

        let mut pixels_to_sample = self.collect_pixels_to_supersample(
            from_row,
            from_col,
            to_row,
            to_col,
            buffer,
            minimum_luminance,
        );

        // Initialize the output buffer with the base 1-spp colors from the raw buffer.
        let mut output_buffer: Vec<Radiance> = buffer.iter().map(|sample| sample.color).collect();

        if let Some(mask_color) = self.scene.sampling_mask_color {
            // The mask is a fixed diagnostic color; pre-divide by the
            // exposure so the later exposure multiply restores it exactly
            // and the overlay looks identical at every --exposure.
            let mask_color = mask_color.mul_color_part(1.0 / self.exposure);
            for pixel in &pixels_to_sample {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = Radiance::from_straight(mask_color);
            }
        } else {
            self.supersample(samples_per_axis, &mut pixels_to_sample)?;
        }

        for pixel in pixels_to_sample {
            if let Some(sample_color) = pixel.result {
                let pixel_index = self.get_pixel_index(
                    pixel.row,
                    pixel.col,
                    to_col - from_col,
                    from_row,
                    from_col,
                );
                output_buffer[pixel_index] = sample_color;
            }
        }

        let width = (to_col - from_col) as usize;
        let height = (to_row - from_row) as usize;
        self.add_star_layer(buffer, width, height, &mut output_buffer)?;

        info!(
            "Finished rendering section from ({}, {}) to ({}, {})",
            from_row, from_col, to_row, to_col
        );
        Ok(output_buffer)
    }

    fn supersample(
        &self,
        samples_per_axis: usize,
        pixels_to_sample: &mut Vec<PixelToSample>,
    ) -> Result<(), RaytracerError> {
        info!("Supersampling {} pixels", pixels_to_sample.len());

        let count = AtomicUsize::new(0);
        use indicatif::{ProgressBar, ProgressStyle};
        let pb = ProgressBar::new(pixels_to_sample.len() as u64);
        pb.set_style(ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
            .map_err(RaytracerError::ProgressBarTemplateError)?
            .progress_chars("█▇▆▅▄▃▂▁  "));

        pixels_to_sample.par_iter_mut().for_each(|pixel| {
            count.fetch_add(1, Ordering::SeqCst);
            pb.set_position(count.load(Ordering::Relaxed) as u64);

            let mut sample_colors = RadianceMean::default();
            for stratum_row in 0..samples_per_axis {
                for stratum_col in 0..samples_per_axis {
                    let (dx, dy) = stratified_sample_offset(
                        pixel.row as i64,
                        pixel.col as i64,
                        stratum_row,
                        stratum_col,
                        samples_per_axis,
                    );
                    let ray = self.scene.camera.get_ray_for_offset(
                        pixel.row as i64,
                        pixel.col as i64,
                        dx,
                        dy,
                    );
                    match self.scene.color_of_ray(&ray) {
                        Ok(sample) => {
                            sample_colors.add(sample.color, 1.0);
                        }
                        Err(err) => {
                            error!(
                                "Unable to compute color for ray at pixel ({}, {}): {:?}",
                                pixel.col, pixel.row, err
                            );
                        }
                    }
                }
            }
            // Divide by the number of rays that actually returned a colour, so
            // a failed sub-sample does not bias the pixel toward black. If all
            // failed, leave result = None and keep the base 1-spp colour.
            pixel.result = sample_colors.mean();
        });
        pb.finish();

        Ok(())
    }

    fn collect_pixels_to_supersample(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        buffer: &[RaySample],
        minimum_luminance: f64,
    ) -> Vec<PixelToSample> {
        let mut pixels_to_sample: Vec<PixelToSample> = Vec::new();
        for row in from_row..to_row {
            for col in from_col..to_col {
                let pixel_index =
                    self.get_pixel_index(row, col, to_col - from_col, from_row, from_col);
                for (row_shift, col_shift) in [
                    (-1, -1),
                    (-1, 0),
                    (-1, 1),
                    (0, -1),
                    (0, 1),
                    (1, -1),
                    (1, 0),
                    (1, 1),
                ] {
                    let neighbor_row = row as i32 + row_shift;
                    let neighbor_col = col as i32 + col_shift;

                    // Neighbours outside this section are skipped. For a full
                    // render that is only the true image border, and for a
                    // standalone section crop (previewing a region) it just
                    // leaves the crop's 1px border compared within the crop,
                    // which is cosmetically negligible. It would only matter if
                    // separately-rendered sections were stitched into one image,
                    // where an edge crossing a seam could stay at 1 spp and show
                    // a line; making that seam-free would need a 1px selection
                    // halo around the section. Not needed for crop previews.
                    if neighbor_row < from_row as i32
                        || neighbor_row >= to_row as i32
                        || neighbor_col < from_col as i32
                        || neighbor_col >= to_col as i32
                    {
                        continue;
                    }
                    let neighbor_index = self.get_pixel_index(
                        neighbor_row as u32,
                        neighbor_col as u32,
                        to_col - from_col,
                        from_row,
                        from_col,
                    );

                    if let (Some(pixel), Some(neighbor)) =
                        (buffer.get(pixel_index), buffer.get(neighbor_index))
                        && should_supersample_pair(
                            pixel,
                            neighbor,
                            &self.scene.adaptive_sampling,
                            minimum_luminance,
                        )
                    {
                        pixels_to_sample.push(PixelToSample {
                            row,
                            col,
                            result: None,
                        });
                        break;
                    }
                }
            }
        }

        pixels_to_sample
    }

    pub fn render_section(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        filename: String,
    ) -> Result<(), RaytracerError> {
        if filename.ends_with(".hdr") {
            info!("Creating HDR image");
            let raw_cie =
                self.render_section_to_radiance_buffer(from_row, from_col, to_row, to_col)?;
            let buffer: Vec<f32> = raw_cie
                .into_iter()
                .map(|c| xyz_to_linear_srgb(&c.to_xyz_over_black()))
                .flat_map(|c| {
                    [
                        c.x.max(0.0) as f32,
                        c.y.max(0.0) as f32,
                        c.z.max(0.0) as f32,
                    ]
                })
                .collect();
            let imgbuf_hdr: ImageBuffer<Rgb<f32>, Vec<f32>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf_hdr
                .save_with_format(&filename, ImageFormat::Hdr)
                .map_err(RaytracerError::ImageError)?;
        } else {
            info!("Creating non-HDR image");
            info!(
                "Tone mapping method: {:?}, exposure: {}",
                self.tone_mapping, self.exposure
            );
            let cie_pixels =
                self.render_section_to_radiance_buffer(from_row, from_col, to_row, to_col)?;
            let cie_pixels = cie_pixels
                .into_iter()
                .map(Radiance::to_xyz_over_black)
                .collect();
            let linear_srgb = xyz_to_linear_srgb_buffer(&cie_pixels);
            let colors = linear_srgb_to_srgb_buffer(&linear_srgb, self.exposure, self.tone_mapping);
            let buffer: Vec<u8> = colors.iter().flat_map(|c| [c.r, c.g, c.b]).collect();
            let imgbuf: ImageBuffer<Rgb<u8>, Vec<u8>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf.save(&filename).map_err(RaytracerError::ImageError)?;
        }

        info!("saved image to {}", filename);
        Ok(())
    }

    pub fn integrate_ray_at_point(
        &self,
        row: i64,
        col: i64,
    ) -> Result<(IntegratedRay, Option<StopReason>), RaytracerError> {
        let ray = self.scene.camera.get_ray_for(row, col);
        info!("ray for {}-{} is: {:?}", row, col, ray);
        self.scene.integrate_ray(&ray)
    }
}

#[cfg(test)]
mod tube_tests;

#[cfg(test)]
mod compositing_tests;

#[cfg(test)]
mod tests {
    use super::{
        MICHELSON_DENOMINATOR_EPSILON, luminance_contrast, should_supersample_pair,
        stratified_sample_offset,
    };
    use crate::configuration::AdaptiveSamplingConfig;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::radiance::Radiance;
    use crate::rendering::ray::Ray;
    use crate::rendering::scene::{EscapeInfo, RayClass, RaySample};

    fn sample(y: f64, alpha: f64, ray_class: RayClass) -> RaySample {
        RaySample {
            color: Radiance::new(0.0, y, 0.0, 1.0 - alpha),
            ray_class,
            accumulated_angular_distance: 0.0,
            ray: Ray::new(
                0,
                0,
                Point::new(
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    crate::geometry::point::CoordinateSystem::Cartesian,
                ),
                FourVector::new_cartesian(0.0, 0.0, 1.0, 0.0),
            ),
        }
    }

    #[test]
    fn stratified_offsets_stay_in_their_cells() {
        let samples_per_axis = 4;
        for stratum_row in 0..samples_per_axis {
            for stratum_col in 0..samples_per_axis {
                let (dx, dy) =
                    stratified_sample_offset(17, 23, stratum_row, stratum_col, samples_per_axis);
                let cell_size = 1.0 / samples_per_axis as f64;
                assert!(
                    (stratum_col as f64 * cell_size..(stratum_col + 1) as f64 * cell_size)
                        .contains(&dx)
                );
                assert!(
                    (stratum_row as f64 * cell_size..(stratum_row + 1) as f64 * cell_size)
                        .contains(&dy)
                );
                assert_eq!(
                    (dx, dy),
                    stratified_sample_offset(17, 23, stratum_row, stratum_col, samples_per_axis,)
                );
            }
        }
    }

    #[test]
    fn michelson_contrast_uses_the_named_epsilon() {
        let black = Radiance::BLACK;
        let faint = Radiance::new(0.0, MICHELSON_DENOMINATOR_EPSILON, 0.0, 0.0);
        assert_eq!(luminance_contrast(&black, &black), 0.0);
        assert_eq!(luminance_contrast(&black, &faint), 0.5);
    }

    #[test]
    fn class_boundaries_are_always_supersampled() {
        let config = AdaptiveSamplingConfig::default();
        let escaped = sample(
            0.0,
            1.0,
            RayClass::Escaped(EscapeInfo {
                x: 0.0,
                y: 0.0,
                z: 1.0,
                redshift: 1.0,
            }),
        );
        let captured = sample(0.0, 1.0, RayClass::Captured);

        assert!(should_supersample_pair(&escaped, &captured, &config, 100.0));
        assert!(should_supersample_pair(&captured, &escaped, &config, 100.0));
    }

    #[test]
    fn background_contrast_does_not_trigger_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.0,
            opacity_contrast_threshold: 0.0,
            ..Default::default()
        };
        let dark = sample(
            1.0,
            0.0,
            RayClass::Escaped(EscapeInfo {
                x: 0.0,
                y: 0.0,
                z: 1.0,
                redshift: 1.0,
            }),
        );
        let bright = sample(
            100.0,
            1.0,
            RayClass::Escaped(EscapeInfo {
                x: 1.0,
                y: 0.0,
                z: 0.0,
                redshift: 1.0,
            }),
        );

        assert!(!should_supersample_pair(&dark, &bright, &config, 0.0));
    }

    #[test]
    fn visible_object_contrast_triggers_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.2,
            opacity_contrast_threshold: 0.2,
            ..Default::default()
        };

        assert!(should_supersample_pair(
            &sample(2.0, 1.0, RayClass::Hit),
            &sample(1.0, 1.0, RayClass::Hit),
            &config,
            1.0,
        ));
        assert!(should_supersample_pair(
            &sample(2.0, 0.6, RayClass::Hit),
            &sample(2.0, 0.9, RayClass::Hit),
            &config,
            1.0,
        ));
    }

    #[test]
    fn faint_object_contrast_does_not_trigger_supersampling() {
        let config = AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.0,
            opacity_contrast_threshold: 0.0,
            ..Default::default()
        };

        assert!(!should_supersample_pair(
            &sample(1.0, 0.0, RayClass::Hit),
            &sample(0.0, 1.0, RayClass::Hit),
            &config,
            1.0,
        ));
    }
}
