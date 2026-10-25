//! Convergence check for a finished render: re-trace and re-shade a sample of
//! pixels at a tighter integrator tolerance and count the ones that change.
//!
//! Two things are compared, because they fail differently. Where the geodesic
//! ends: near the photon orbit a loose tolerance changes how many times a ray
//! loops around the hole and sends it to another patch of sky. And the shaded
//! colour (objects and celestial texture, not catalogue stars): a path that
//! ends in the same place can still cross a volumetric disc somewhere else,
//! and with fine noise textures that changes the pixel. On the README hero at
//! epsilon 1e-5 the second is by far the larger effect.

use super::Raytracer;
use crate::geometry::geometry::Geometry;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::StopReason;
use crate::rendering::ray::IntegratedRay;
use log::{info, warn};
use rayon::prelude::*;
use serde::Serialize;

/// How much tighter the reference tolerance is than the render's.
pub const REFERENCE_FACTOR: f64 = 100.0;

/// The reference never goes below this: the error control is absolute, and
/// with coordinates up to max_radius (~1e4) a per-step error much under 1e-12
/// is below double precision, so the step size collapses and the reference
/// runs out of steps instead of getting more accurate.
pub const MIN_REFERENCE_EPSILON: f64 = 1e-12;

/// The reference gets this many times the render's step budget.
pub const REFERENCE_STEP_FACTOR: usize = 100;

/// Share of rays ending differently above which the check warns.
pub const WARN_FRACTION: f64 = 0.01;

/// Share of lit samples changing colour above which the check warns. Higher
/// than WARN_FRACTION because single unsupersampled rays through a volumetric
/// disc's fine noise flicker more than the final image does: on the README
/// hero (KerrBL, 320x180, 1000 samples) the 1e-9 default sits at 1.9%, 1e-7 at
/// 7.2% and 1e-5 at 18.4%.
pub const COLOUR_WARN_FRACTION: f64 = 0.05;

/// Relative luminance difference at which a lit sample counts as changed.
pub const COLOUR_TOLERANCE: f64 = 0.25;

/// Samples dimmer than this fraction of the brightest reference sample are
/// left out of the colour comparison (black sky, deep shadow).
const LIT_FLOOR: f64 = 1e-3;

#[derive(Debug, Clone, Serialize)]
pub struct AccuracyReport {
    pub samples: usize,
    pub epsilon: f64,
    pub reference_epsilon: f64,
    /// Escape directions further apart than this count as changed: about one
    /// pixel at the image centre.
    pub threshold_rad: f64,
    /// Pixels whose ray stops for a different reason (e.g. escapes instead of
    /// falling in) or escapes more than `threshold_rad` away.
    pub changed: usize,
    pub changed_fraction: f64,
    /// Largest escape-direction difference among rays that escape both times.
    pub max_deviation_rad: f64,
    /// Samples bright enough to compare colours (see LIT_FLOOR).
    pub lit_samples: usize,
    /// Lit samples whose luminance differs by more than COLOUR_TOLERANCE.
    pub colour_changed: usize,
    pub colour_changed_fraction: f64,
}

/// Where a traced ray ends up: its stop reason and, if it reached the
/// celestial sphere, the unit direction it escaped in.
fn outcome<C: crate::geometry::coordinates::conversion::CartesianConvertible>(
    ray: &IntegratedRay<C>,
    stop: Option<StopReason>,
) -> (Option<StopReason>, Option<nalgebra::Vector3<f64>>) {
    let direction = match stop {
        Some(StopReason::CelestialSphereReached) => ray
            .steps
            .last()
            .map(|s| s.x.get_spatial_vector_cartesian())
            .filter(|v| v.norm() > 0.0)
            .map(|v| v.normalize()),
        _ => None,
    };
    (stop, direction)
}

/// `samples` pixel positions spread evenly over the section, row-major.
fn sample_grid(
    from_row: u32,
    from_col: u32,
    to_row: u32,
    to_col: u32,
    samples: usize,
) -> Vec<(i64, i64)> {
    let rows = (to_row - from_row) as f64;
    let cols = (to_col - from_col) as f64;
    if samples == 0 || rows <= 0.0 || cols <= 0.0 {
        return Vec::new();
    }
    // An n_r x n_c grid with roughly the section's aspect ratio.
    let n_c = ((samples as f64 * cols / rows).sqrt().round() as usize).max(1);
    let n_r = samples.div_ceil(n_c).max(1);
    let mut grid = Vec::with_capacity(n_r * n_c);
    for i in 0..n_r {
        for j in 0..n_c {
            let row = from_row as f64 + (i as f64 + 0.5) * rows / n_r as f64;
            let col = from_col as f64 + (j as f64 + 0.5) * cols / n_c as f64;
            grid.push((row as i64, col as i64));
        }
    }
    grid.truncate(samples);
    grid
}

/// Half the samples on an even grid, half spread over the lit pixels of the
/// rendered section (`luminance` is row-major over it). A frame is often mostly
/// black sky, where nothing can change colour, so a grid alone puts few
/// samples on the disc.
fn sample_pixels(
    from_row: u32,
    from_col: u32,
    to_row: u32,
    to_col: u32,
    luminance: &[f64],
    samples: usize,
) -> Vec<(i64, i64)> {
    let width = (to_col - from_col) as usize;
    let mut pixels = sample_grid(from_row, from_col, to_row, to_col, samples / 2);
    let brightest = luminance.iter().copied().fold(0.0, f64::max);
    let lit: Vec<usize> = (0..luminance.len())
        .filter(|&i| luminance[i] > LIT_FLOOR * brightest)
        .collect();
    let wanted = samples - pixels.len();
    if width > 0 && !lit.is_empty() && wanted > 0 {
        let stride = lit.len() as f64 / wanted as f64;
        for k in 0..wanted.min(lit.len()) {
            let i = lit[((k as f64 + 0.5) * stride) as usize];
            pixels.push((
                from_row as i64 + (i / width) as i64,
                from_col as i64 + (i % width) as i64,
            ));
        }
    }
    pixels.sort_unstable();
    pixels.dedup();
    pixels
}

impl<G: Geometry> Raytracer<'_, G> {
    /// Re-trace and re-shade about `samples` pixels of the section at
    /// `epsilon / REFERENCE_FACTOR` (see `sample_pixels` for which) and compare
    /// where their geodesics end and how bright they are. Logs the result and
    /// warns when more than `WARN_FRACTION` of either changed.
    pub fn check_accuracy(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        luminance: &[f64],
        samples: usize,
    ) -> Result<AccuracyReport, RaytracerError> {
        let integrator = &self.scene.integrator;
        let epsilon = integrator.epsilon();
        let reference_epsilon = (epsilon / REFERENCE_FACTOR).max(MIN_REFERENCE_EPSILON);
        let reference = integrator.with_accuracy(
            reference_epsilon,
            integrator.max_steps().saturating_mul(REFERENCE_STEP_FACTOR),
        );
        let threshold_rad = self.scene.camera.pixel_angle();
        let pixels = sample_pixels(from_row, from_col, to_row, to_col, luminance, samples);

        info!(
            "Accuracy check: re-tracing {} pixels at epsilon {:e} (render: {:e})",
            pixels.len(),
            reference_epsilon,
            epsilon
        );
        // (geodesic changed, escape deviation, luminance, reference luminance)
        let results: Vec<(bool, Option<f64>, f64, f64)> = pixels
            .par_iter()
            .map(|&(row, col)| {
                let ray = self.scene.camera.get_ray_for(row, col);
                let (a, stop_a) = integrator.integrate(&ray)?;
                let (b, stop_b) = reference.integrate(&ray)?;
                let (reason_a, dir_a) = outcome(&a, stop_a);
                let (reason_b, dir_b) = outcome(&b, stop_b);
                let deviation = match (dir_a, dir_b) {
                    (Some(u), Some(v)) => Some(u.dot(&v).clamp(-1.0, 1.0).acos()),
                    _ => None,
                };
                let changed = reason_a != reason_b || deviation.is_some_and(|d| d > threshold_rad);
                let y = self.scene.color_of_ray_with(integrator, &ray)?.color.y;
                let y_ref = self.scene.color_of_ray_with(&reference, &ray)?.color.y;
                Ok((changed, deviation, y, y_ref))
            })
            .collect::<Result<_, RaytracerError>>()?;

        let changed = results.iter().filter(|r| r.0).count();
        let max_deviation_rad = results.iter().filter_map(|r| r.1).fold(0.0, f64::max);
        let brightest = results.iter().map(|r| r.3).fold(0.0, f64::max);
        let lit: Vec<_> = results
            .iter()
            .filter(|r| r.3 > LIT_FLOOR * brightest)
            .collect();
        let colour_changed = lit
            .iter()
            .filter(|r| (r.2 - r.3).abs() > COLOUR_TOLERANCE * r.3)
            .count();
        let report = AccuracyReport {
            samples: results.len(),
            epsilon,
            reference_epsilon,
            threshold_rad,
            changed,
            changed_fraction: if results.is_empty() {
                0.0
            } else {
                changed as f64 / results.len() as f64
            },
            max_deviation_rad,
            lit_samples: lit.len(),
            colour_changed,
            colour_changed_fraction: if lit.is_empty() {
                0.0
            } else {
                colour_changed as f64 / lit.len() as f64
            },
        };

        let summary = format!(
            "Accuracy check at epsilon {:e}: {} of {} sampled rays ({:.2}%) end differently \
             (threshold {:.2e} rad, max escape deviation {:.2e} rad); {} of {} lit samples \
             ({:.2}%) change luminance by more than {:.0}%",
            reference_epsilon,
            report.changed,
            report.samples,
            100.0 * report.changed_fraction,
            threshold_rad,
            max_deviation_rad,
            report.colour_changed,
            report.lit_samples,
            100.0 * report.colour_changed_fraction,
            100.0 * COLOUR_TOLERANCE
        );
        if report.changed_fraction > WARN_FRACTION
            || report.colour_changed_fraction > COLOUR_WARN_FRACTION
        {
            warn!("{summary}; the render is not converged, lower --epsilon");
        } else {
            info!("{summary}");
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::{sample_grid, sample_pixels};

    #[test]
    fn sample_pixels_puts_half_on_lit_pixels() {
        // A 10x10 section, dark except for a bright 2x2 block at rows 4-5,
        // cols 6-7 (section-relative).
        let mut luminance = vec![0.0; 100];
        for (r, c) in [(4, 6), (4, 7), (5, 6), (5, 7)] {
            luminance[r * 10 + c] = 1.0;
        }
        let pixels = sample_pixels(20, 30, 30, 40, &luminance, 8);
        let on_block = pixels
            .iter()
            .filter(|&&(r, c)| (24..26).contains(&r) && (36..38).contains(&c))
            .count();
        assert_eq!(on_block, 4, "{pixels:?}");
        assert!(
            pixels
                .iter()
                .all(|&(r, c)| (20..30).contains(&r) && (30..40).contains(&c))
        );
    }

    #[test]
    fn sample_grid_spreads_over_the_section() {
        let grid = sample_grid(100, 200, 300, 600, 32);
        assert_eq!(grid.len(), 32);
        assert!(
            grid.iter()
                .all(|&(r, c)| (100..300).contains(&r) && (200..600).contains(&c))
        );
        // Wider than tall, so more columns than rows in the grid.
        let distinct_cols = grid
            .iter()
            .map(|&(_, c)| c)
            .collect::<std::collections::BTreeSet<_>>();
        let distinct_rows = grid
            .iter()
            .map(|&(r, _)| r)
            .collect::<std::collections::BTreeSet<_>>();
        assert!(distinct_cols.len() > distinct_rows.len());
    }

    #[test]
    fn sample_grid_handles_empty_requests() {
        assert!(sample_grid(0, 0, 10, 10, 0).is_empty());
        assert!(sample_grid(5, 5, 5, 10, 4).is_empty());
    }
}
