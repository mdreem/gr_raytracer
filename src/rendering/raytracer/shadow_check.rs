//! Exact check of a finished render's shadow: compare where each pixel's
//! primary ray ended (horizon or celestial sphere) with the geometry's exact
//! classification of the same camera ray (`Geometry::analytic_capture`).
//!
//! Unlike `check_accuracy` this is not self-convergence on a sample: the
//! reference is exact and covers every pixel, and it costs one small
//! eigenvalue problem per pixel. Pixels next to the critical curve are left
//! out, since there the true outcome changes within the pixel: a pixel is
//! compared only when all of its 3x3 neighbours have the same exact class.
//! It checks the geodesics only, not what they cross on the way (disc,
//! catalogue stars).

use super::Raytracer;
use crate::geometry::geometry::Geometry;
use crate::rendering::integration::integrator::StopReason;
use log::{info, warn};
use rayon::prelude::*;
use serde::Serialize;

/// How many disagreeing pixels the report lists.
const MAX_LISTED: usize = 10;

#[derive(Debug, Clone, Serialize)]
pub struct ShadowCheckReport {
    /// Pixels whose primary ray the geometry classifies exactly.
    pub classified: usize,
    /// Classified pixels whose 3x3 neighbourhood has the same exact class.
    pub compared: usize,
    /// Compared pixels inside the exact shadow.
    pub inside_shadow: usize,
    /// Compared pixels whose rendered ray ends differently.
    pub disagreements: usize,
    pub disagreement_fraction: f64,
    /// Disagreements that escape although the exact ray falls in.
    pub escaped_inside_shadow: usize,
    /// Disagreements that do not escape although the exact ray does. Includes
    /// rays stopped as NaN, closed orbit or out of steps: the image shows
    /// those as shadow too.
    pub captured_outside_shadow: usize,
    /// Pixels (Chebyshev distance) from the deepest disagreement to the
    /// nearest pixel left out of the comparison, i.e. how far from the
    /// critical curve the render goes wrong.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_depth_px: Option<usize>,
    /// Up to MAX_LISTED disagreeing pixels as [row, col], deepest first.
    pub worst_pixels: Vec<[i64; 2]>,
}

/// Chebyshev distance from each pixel to the nearest pixel with `seed` set
/// (0 on seeds, None if there are none), by the two-pass chamfer transform.
fn distance_to_seeds(seed: &[bool], width: usize, height: usize) -> Vec<Option<usize>> {
    const FAR: usize = usize::MAX / 2;
    let mut d: Vec<usize> = seed.iter().map(|&s| if s { 0 } else { FAR }).collect();
    for row in 0..height {
        for col in 0..width {
            let i = row * width + col;
            for (dr, dc) in [(-1, -1), (-1, 0), (-1, 1), (0, -1)] {
                if let Some(j) = neighbour(row, col, dr, dc, width, height) {
                    d[i] = d[i].min(d[j] + 1);
                }
            }
        }
    }
    for row in (0..height).rev() {
        for col in (0..width).rev() {
            let i = row * width + col;
            for (dr, dc) in [(1, 1), (1, 0), (1, -1), (0, 1)] {
                if let Some(j) = neighbour(row, col, dr, dc, width, height) {
                    d[i] = d[i].min(d[j] + 1);
                }
            }
        }
    }
    d.into_iter().map(|v| (v < FAR).then_some(v)).collect()
}

fn neighbour(
    row: usize,
    col: usize,
    dr: isize,
    dc: isize,
    width: usize,
    height: usize,
) -> Option<usize> {
    let r = row.checked_add_signed(dr).filter(|&r| r < height)?;
    let c = col.checked_add_signed(dc).filter(|&c| c < width)?;
    Some(r * width + c)
}

/// Compare the exact class of each pixel (`Some(true)` = captured) with the
/// rendered one (`true` = did not escape). Both are row-major over a
/// `width` x `height` section starting at `(from_row, from_col)`.
fn compare(
    exact: &[Option<bool>],
    rendered_captured: &[bool],
    width: usize,
    height: usize,
    from_row: u32,
    from_col: u32,
) -> ShadowCheckReport {
    let compared: Vec<bool> = (0..exact.len())
        .map(|i| {
            let (row, col) = (i / width, i % width);
            exact[i].is_some()
                && (-1..=1).all(|dr| {
                    (-1..=1).all(|dc| {
                        neighbour(row, col, dr, dc, width, height)
                            .is_none_or(|j| exact[j] == exact[i])
                    })
                })
        })
        .collect();
    let left_out: Vec<bool> = compared.iter().map(|c| !c).collect();
    let depth = distance_to_seeds(&left_out, width, height);

    let mut disagreeing: Vec<(Option<usize>, usize)> = Vec::new();
    let (mut escaped_inside_shadow, mut captured_outside_shadow) = (0, 0);
    for i in (0..exact.len()).filter(|&i| compared[i]) {
        let captured = exact[i] == Some(true);
        if captured != rendered_captured[i] {
            if captured {
                escaped_inside_shadow += 1;
            } else {
                captured_outside_shadow += 1;
            }
            disagreeing.push((depth[i], i));
        }
    }
    // Deepest first; None (no critical curve in the section) is deepest of all.
    disagreeing.sort_by_key(|&(d, i)| (std::cmp::Reverse(d.unwrap_or(usize::MAX)), i));

    let compared_count = compared.iter().filter(|&&c| c).count();
    ShadowCheckReport {
        classified: exact.iter().filter(|e| e.is_some()).count(),
        compared: compared_count,
        inside_shadow: (0..exact.len())
            .filter(|&i| compared[i] && exact[i] == Some(true))
            .count(),
        disagreements: disagreeing.len(),
        disagreement_fraction: if compared_count == 0 {
            0.0
        } else {
            disagreeing.len() as f64 / compared_count as f64
        },
        escaped_inside_shadow,
        captured_outside_shadow,
        max_depth_px: disagreeing.first().and_then(|&(d, _)| d),
        worst_pixels: disagreeing
            .iter()
            .take(MAX_LISTED)
            .map(|&(_, i)| {
                [
                    from_row as i64 + (i / width) as i64,
                    from_col as i64 + (i % width) as i64,
                ]
            })
            .collect(),
    }
}

impl<G: Geometry> Raytracer<'_, G> {
    /// Classify every pixel's camera ray exactly and compare with
    /// `stop_reasons` (row-major over the section, from `render_section`).
    /// Returns None, after logging why, if the geometry classifies none of
    /// them (no exact solution, or no horizon). Warns on any disagreement.
    pub fn check_shadow(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        stop_reasons: &[Option<StopReason>],
    ) -> Option<ShadowCheckReport> {
        let width = (to_col - from_col) as usize;
        let height = (to_row - from_row) as usize;
        let exact: Vec<Option<bool>> = (0..width * height)
            .into_par_iter()
            .map(|i| {
                let row = from_row as i64 + (i / width) as i64;
                let col = from_col as i64 + (i % width) as i64;
                let ray = self.scene.camera.get_ray_for(row, col);
                self.scene
                    .geometry
                    .analytic_capture(&ray.position, &ray.momentum)
            })
            .collect();
        if exact.iter().all(Option::is_none) {
            info!("Shadow check skipped: no exact capture classification for this geometry");
            return None;
        }
        let rendered_captured: Vec<bool> = stop_reasons
            .iter()
            .map(|&s| s != Some(StopReason::CelestialSphereReached))
            .collect();
        let report = compare(
            &exact,
            &rendered_captured,
            width,
            height,
            from_row,
            from_col,
        );

        let summary = format!(
            "Shadow check: {} of {} pixels compared with the exact shadow ({} inside it), \
             {} disagree ({} escape inside the shadow, {} do not escape outside it)",
            report.compared,
            width * height,
            report.inside_shadow,
            report.disagreements,
            report.escaped_inside_shadow,
            report.captured_outside_shadow,
        );
        if report.disagreements > 0 {
            warn!(
                "{summary}; up to {} px from the critical curve, e.g. [row, col] {:?}; \
                 lower --epsilon",
                report
                    .max_depth_px
                    .map_or("an unknown number of".to_string(), |d| d.to_string()),
                report.worst_pixels,
            );
        } else {
            info!("{summary}");
        }
        Some(report)
    }
}

#[cfg(test)]
mod tests {
    use super::{compare, distance_to_seeds};

    #[test]
    fn distance_is_chebyshev() {
        let mut seed = vec![false; 25];
        seed[0] = true;
        let d = distance_to_seeds(&seed, 5, 5);
        assert_eq!(d[0], Some(0));
        assert_eq!(d[6], Some(1));
        assert_eq!(d[24], Some(4));
        assert_eq!(d[4], Some(4));
        assert!(
            distance_to_seeds(&[false; 4], 2, 2)
                .iter()
                .all(Option::is_none)
        );
    }

    #[test]
    fn compare_skips_the_critical_curve_and_counts_the_rest() {
        // A 6x8 section, captured in columns 0..4, escaping in 4..8.
        let (width, height) = (8, 6);
        let exact: Vec<Option<bool>> = (0..width * height).map(|i| Some(i % width < 4)).collect();
        let mut rendered: Vec<bool> = exact.iter().map(|e| e.unwrap()).collect();
        // Wrong right at the edge (columns 3 and 4): left out, not counted.
        rendered[2 * width + 3] = false;
        rendered[2 * width + 4] = true;
        let report = compare(&exact, &rendered, width, height, 0, 0);
        assert_eq!(report.classified, 48);
        assert_eq!(report.compared, 36);
        assert_eq!(report.inside_shadow, 18);
        assert_eq!(report.disagreements, 0);

        // Escapes two pixels inside the shadow and falls in one outside it.
        rendered[3 * width + 1] = false;
        rendered[4 * width + 6] = true;
        let report = compare(&exact, &rendered, width, height, 100, 200);
        assert_eq!(report.disagreements, 2);
        assert_eq!(report.escaped_inside_shadow, 1);
        assert_eq!(report.captured_outside_shadow, 1);
        assert_eq!(report.max_depth_px, Some(2));
        assert_eq!(report.worst_pixels, vec![[103, 201], [104, 206]]);

        // Unclassified pixels are left out like the curve.
        let mut exact = exact;
        exact[3 * width + 1] = None;
        let report = compare(&exact, &rendered, width, height, 0, 0);
        assert_eq!(report.classified, 47);
        assert_eq!(report.disagreements, 1);
    }
}
