//! Heavy, ignored end-to-end test: a single bright star placed directly behind
//! a Schwarzschild hole must lens into an Einstein ring. This is the in-repo
//! reproduction of the single-star ring figure from the star-catalogue work: a
//! correctness check on the lensing + star gather, and a way to eyeball the
//! ring without the external one-star parquet + render pipeline.
//!
//! Ignored because it runs a real (small) render. Run it with:
//!   cargo test --release ring -- --ignored --nocapture
//! It also writes `target/single_star_ring.pgm` for visual inspection.
//!
//! Set CURVED_MEMBERSHIP=1 to exercise the curved-boundary gather, and
//! POLE_ROT_DEG=<deg> to steer the coordinate pole off the ring.

use super::*;
use crate::geometry::geometry::{HasCoordinateSystem, SupportQuantities};
use crate::geometry::point::Point;
use crate::geometry::schwarzschild::Schwarzschild;
use crate::rendering::camera::Camera;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::octree::Octree;
use crate::rendering::star_catalog::{Star, StarCatalog};
use crate::rendering::texture::{TemperatureData, TextureData, TextureMap, UVCoordinates};
use crate::scene_objects::objects::Objects;
use nalgebra::Vector3;
use std::f64::consts::PI;
use std::sync::Arc;

/// A pure-black background so only the lensed star contributes.
struct BlackSky;

impl TextureMap for BlackSky {
    fn color_at_uv(
        &self,
        _: &UVCoordinates,
        _: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
        Ok(CIETristimulus::new(0.0, 0.0, 0.0, 1.0))
    }
}

/// One bright star on the sky direction directly behind a camera on the -y
/// axis, i.e. `(0, 1, 0)` (ra = 90 deg, dec = 0). Perfect alignment, so its
/// image is a full Einstein ring.
fn star_directly_behind() -> Star {
    Star {
        source_id: 1,
        ra_deg: 90.0,
        dec_deg: 0.0,
        g_mag: 0.0,
        bp_mag: 0.0,
        rp_mag: 0.0,
        bp_rp: 0.0,
        direction: Vector3::new(0.0, 1.0, 0.0),
        temperature: 5772.0,
        // Y carries the flux; hue is irrelevant for the ring geometry.
        emission_xyz: CIETristimulus::new(0.95, 1.0, 1.09, 1.0),
    }
}

/// Azimuthally-averaged luminance profile about the image centre.
fn radial_profile(buffer: &[Radiance], width: usize, height: usize) -> Vec<f64> {
    let (cx, cy) = ((width as f64 - 1.0) / 2.0, (height as f64 - 1.0) / 2.0);
    let rmax = width.min(height) / 2;
    let mut sum = vec![0.0f64; rmax];
    let mut count = vec![0u32; rmax];
    for row in 0..height {
        for col in 0..width {
            let dr = (((col as f64) - cx).powi(2) + ((row as f64) - cy).powi(2)).sqrt();
            let k = dr as usize;
            if k < rmax {
                sum[k] += buffer[row * width + col].y;
                count[k] += 1;
            }
        }
    }
    (0..rmax)
        .map(|k| if count[k] > 0 { sum[k] / count[k] as f64 } else { 0.0 })
        .collect()
}

/// Write a gamma-corrected grayscale PGM for eyeballing the ring.
fn write_pgm(buffer: &[Radiance], width: usize, height: usize, path: &str) {
    let max = buffer.iter().map(|r| r.y).fold(0.0f64, f64::max).max(1e-12);
    let mut out = format!("P5\n{width} {height}\n255\n").into_bytes();
    for pixel in buffer.iter() {
        let v = (pixel.y / max).clamp(0.0, 1.0).powf(1.0 / 2.2);
        out.push((v * 255.0 + 0.5) as u8);
    }
    let _ = std::fs::write(path, out);
}

#[test]
#[ignore = "heavy render; run with `cargo test --release ring -- --ignored`"]
fn single_star_directly_behind_forms_einstein_ring() {
    let (width, height) = (100usize, 100usize);
    let geometry = Schwarzschild::new(1.0, 1e-4);

    // Honour POLE_ROT_DEG before any coordinate conversion below: the first
    // cartesian->spherical conversion locks the frame-rotation OnceLock, so
    // setting it afterwards would be silently ignored.
    if let Some(deg) = std::env::var("POLE_ROT_DEG").ok().and_then(|s| s.parse().ok()) {
        crate::geometry::spherical_coordinates_helper::set_frame_rotation_deg(deg);
    }

    // Camera on the -y axis, looking back at the hole (theta = pi), so the
    // shadow and ring are centred in the frame. Use the geometry's normalized
    // static-observer 4-velocity, not a raw (1,0,0,0), so the tetrad is valid.
    let position = Point::new_cartesian(0.0, 0.0, -20.0, 0.0)
        .to_coordinate_system(geometry.coordinate_system());
    let velocity = geometry.get_stationary_velocity_at(&position);
    let camera = Camera::new(
        position,
        velocity,
        PI / 6.0,
        height as i64,
        width as i64,
        0.0,
        PI,
        0.0,
        &geometry,
    )
    .expect("camera");

    let scene = Scene::new(
        IntegrationConfiguration::new(12000, 400.0, 0.02, 1e-6),
        Objects::new(&geometry),
        TextureData {
            celestial_map: Arc::new(BlackSky),
        },
        &geometry,
        camera,
        false,
        0.0,
    );
    let mut renderer = Raytracer::new(scene, ToneMappingMethod::Reinhard, 1.0);
    // Take the direct star-layer path (no supersampling) so the gather runs.
    renderer.scene.adaptive_sampling.enabled = false;
    renderer.scene.star_catalog = Some(StarCatalog {
        stars: Octree::new(vec![star_directly_behind()]),
    });
    // The production toggles are the --curved-star-membership and
    // --pole-rotation-deg flags; honour the same knobs here via env vars so
    // this test can exercise both paths without a CLI. (POLE_ROT_DEG is applied
    // above, before the coordinate conversion.)
    renderer.scene.curved_star_membership = std::env::var_os("CURVED_MEMBERSHIP").is_some();

    let buffer = renderer
        .render_section_to_radiance_buffer(0, 0, height as u32, width as u32)
        .expect("render");

    write_pgm(&buffer, width, height, "target/single_star_ring.pgm");

    let profile = radial_profile(&buffer, width, height);
    let (peak_radius, &peak) = profile
        .iter()
        .enumerate()
        .skip(4) // ignore the shadow core
        .max_by(|a, b| a.1.total_cmp(b.1))
        .expect("profile");
    let core = profile[..4].iter().copied().fold(0.0, f64::max);

    // The brightness is a ring at a nonzero radius, not a central blob, and the
    // shadow core is dark compared with the ring.
    assert!(
        peak_radius > 6,
        "brightness peak should be an off-centre ring, got radius {peak_radius}"
    );
    assert!(
        peak > 10.0 * core.max(1e-9),
        "ring ({peak:.3e}) should be far brighter than the shadow core ({core:.3e})"
    );
}
