use super::*;
use crate::geometry::euclidean::EuclideanSpace;
use crate::rendering::camera::Camera;
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::star_catalog::StarCatalog;
use crate::rendering::texture::{TemperatureData, TextureData, TextureMap, UVCoordinates};
use crate::rendering::tubetracer::TubeScreenBounds;
use crate::scene_objects::objects::Objects;
use approx::assert_abs_diff_eq;
use std::sync::Arc;

struct CountingTexture {
    calls: Arc<AtomicUsize>,
    fail: bool,
}

impl TextureMap for CountingTexture {
    fn color_at_uv(
        &self,
        _: &UVCoordinates,
        _: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        if self.fail {
            Err(RaytracerError::BelowRISCO)
        } else {
            Ok(CIETristimulus::new(2.0, 1.0, 0.5, 1.0))
        }
    }
}

fn renderer<'a>(
    g: &'a EuclideanSpace,
    calls: &Arc<AtomicUsize>,
    fail: bool,
) -> Raytracer<'a, EuclideanSpace> {
    let camera = Camera::new(
        Point::new_cartesian(0.0, 0.0, 0.0, 10.0),
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
        0.7,
        32,
        32,
        0.0,
        0.0,
        0.0,
        g,
    )
    .unwrap();
    let scene = Scene::new(
        IntegrationConfiguration::new(100, 20.0, 0.01, 1e-8),
        Objects::new(g),
        TextureData {
            celestial_map: Arc::new(CountingTexture {
                calls: calls.clone(),
                fail,
            }),
        },
        g,
        camera,
        false,
        0.0,
    );
    Raytracer::new(scene, ToneMappingMethod::Reinhard, 1.0)
}

fn corners(renderer: &Raytracer<EuclideanSpace>, row: i64, col: i64) -> [RaySample; 4] {
    [
        (row, col),
        (row, col + 1),
        (row + 1, col),
        (row + 1, col + 1),
    ]
    .map(|(r, c)| {
        renderer
            .scene
            .color_of_ray(&renderer.scene.camera.get_ray_for(r, c))
            .unwrap()
    })
}

fn tube(samples: &[RaySample; 4]) -> SampleTube<'_> {
    SampleTube::from_buffer(samples, 0, 0, 2).unwrap()
}

fn sky_area(tube: &SampleTube) -> f64 {
    let escaped =
        [tube.a, tube.b, tube.c, tube.d].map(|sample| sample.ray_class.escaped().unwrap());
    compute_traced_tube_solid_angle(&escaped[0], &escaped[1], &escaped[2], &escaped[3]).unwrap()
}

#[test]
fn tube_subdivision_traces_fractional_positions_and_shares_boundaries() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let renderer = renderer(&g, &calls, false);
    // Absolute positions also exercise a cropped buffer's screen origin.
    let parent_samples = corners(&renderer, 10, 20);
    let parent = tube(&parent_samples);
    calls.store(0, Ordering::Relaxed);
    let new_samples = renderer.trace_subdivision(&parent).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), 5);

    // Offsets relative to pixel (10,20): center, top, left, right, bottom.
    // Right and bottom must use 1.5, not the center's 1.0.
    for (sample, (dx, dy)) in
        new_samples
            .iter()
            .zip([(1.0, 1.0), (1.0, 0.5), (0.5, 1.0), (1.5, 1.0), (1.0, 1.5)])
    {
        let expected = renderer.scene.camera.get_ray_for_offset(10, 20, dx, dy);
        assert_abs_diff_eq!(
            sample.ray.momentum.vector,
            expected.momentum.vector,
            epsilon = 1e-14
        );
        assert_eq!((sample.ray.row, sample.ray.col), (10, 20));
    }

    let children = parent.children(&new_samples);
    assert_eq!(
        children.map(|child| child.screen_bounds),
        [
            TubeScreenBounds {
                top: 10.0,
                bottom: 10.5,
                left: 20.0,
                right: 20.5
            },
            TubeScreenBounds {
                top: 10.0,
                bottom: 10.5,
                left: 20.5,
                right: 21.0
            },
            TubeScreenBounds {
                top: 10.5,
                bottom: 11.0,
                left: 20.0,
                right: 20.5
            },
            TubeScreenBounds {
                top: 10.5,
                bottom: 11.0,
                left: 20.5,
                right: 21.0
            },
        ]
    );
    let children = parent.children(&new_samples);
    assert!(std::ptr::eq(children[0].b, children[1].a));
    assert!(std::ptr::eq(children[0].c, children[2].a));
    assert!(std::ptr::eq(children[0].d, children[3].a));
    assert!(std::ptr::eq(children[1].d, children[3].b));
    assert!(std::ptr::eq(children[2].d, children[3].c));
    assert!(children.iter().all(|child| sky_area(child) > 0.0));
    let screen_area = |tube: &SampleTube| {
        (tube.screen_bounds.bottom - tube.screen_bounds.top)
            * (tube.screen_bounds.right - tube.screen_bounds.left)
    };
    // Screen rectangles partition exactly. Their spherical chord areas
    // only approximate the camera's curved stereographic boundaries.
    assert_eq!(
        children.iter().map(screen_area).sum::<f64>(),
        screen_area(&parent)
    );

    // Integer pixel labels are identical, but a second split must still
    // refine to the lower-right child's own fractional center.
    let grandchildren = renderer.trace_subdivision(&children[3]).unwrap();
    let expected = renderer.scene.camera.get_ray_for_offset(10, 20, 1.25, 1.25);
    assert_abs_diff_eq!(
        grandchildren[0].ray.momentum.vector,
        expected.momentum.vector,
        epsilon = 1e-14
    );
    assert_eq!(calls.load(Ordering::Relaxed), 10);
}

#[test]
fn tube_subdivision_preserves_flat_space_star_flux() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    let mut samples = corners(&renderer, 10, 20);
    let ray = renderer.scene.camera.get_ray_for_offset(10, 20, 0.67, 0.73);
    let star = Star {
        source_id: 1,
        ra_deg: 0.0,
        dec_deg: 0.0,
        g_mag: 0.0,
        bp_mag: 0.0,
        rp_mag: 0.0,
        bp_rp: 0.0,
        direction: ray.momentum.get_cartesian_vector(&ray.position).normalize(),
        temperature: 5772.0,
        emission_xyz: CIETristimulus::new(2.0, 1.0, 0.5, 1.0),
    };
    renderer.scene.star_catalog = Some(StarCatalog { stars: vec![star] });
    // Force the real winding-triggered branch. One corner retains the
    // winding discontinuity into successive children, exercising depth.
    samples[0].accumulated_angular_distance = 10.0;
    renderer.scene.max_subdivision_depth = 0;
    let expected = renderer.handle_tube(&tube(&samples)).unwrap();
    assert_abs_diff_eq!(expected.y, 1.0, epsilon = 1e-12);
    for depth in [1, 2, 3] {
        renderer.scene.max_subdivision_depth = depth;
        let actual = renderer.handle_tube(&tube(&samples)).unwrap();
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-10);
        assert!(actual.iter().all(|v| v.is_finite()));
    }
}

#[test]
fn tube_subdivision_depth_limit_counts_actual_splits() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    let samples = corners(&renderer, 0, 0);
    let parent = tube(&samples);
    calls.store(0, Ordering::Relaxed);
    renderer.scene.max_subdivision_depth = 0;
    assert!(renderer.subdivide(&parent, 0).unwrap().is_none());
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    renderer.scene.max_subdivision_depth = 1;
    assert!(renderer.subdivide(&parent, 0).unwrap().is_some());
    assert_eq!(calls.load(Ordering::Relaxed), 5);
    assert!(renderer.subdivide(&parent, 1).unwrap().is_none());
    assert_eq!(calls.load(Ordering::Relaxed), 5);
}

#[test]
fn star_magnification_and_foreground_transmittance_are_applied_once() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    let mut samples = corners(&renderer, 10, 20);
    let star_ray = renderer.scene.camera.get_ray_for_offset(10, 20, 0.67, 0.73);
    renderer.scene.star_catalog = Some(StarCatalog {
        stars: vec![Star {
            source_id: 1,
            ra_deg: 0.0,
            dec_deg: 0.0,
            g_mag: 0.0,
            bp_mag: 0.0,
            rp_mag: 0.0,
            bp_rp: 0.0,
            direction: star_ray
                .momentum
                .get_cartesian_vector(&star_ray.position)
                .normalize(),
            temperature: 5772.0,
            emission_xyz: CIETristimulus::new(2.0, 1.0, 0.5, 1.0),
        }],
    });
    renderer.scene.max_subdivision_depth = 0;
    // A smaller image-side footprint with unchanged escaped corners gives
    // magnification near 1/4. Before C02 that ratio also became star alpha,
    // causing compositing to multiply the demagnified flux a second time.
    for (sample, (dx, dy)) in
        samples
            .iter_mut()
            .zip([(0.75, 0.75), (1.25, 0.75), (0.75, 1.25), (1.25, 1.25)])
    {
        sample.ray.momentum = renderer
            .scene
            .camera
            .get_ray_for_offset(10, 20, dx, dy)
            .momentum;
    }
    let flux = renderer.handle_tube(&tube(&samples)).unwrap();
    assert!(flux.y > 0.24 && flux.y < 0.26);
    for transmittance in [0.0, 0.1, 0.4, 0.75, 1.0] {
        let foreground = Radiance::new(0.2, 0.1, 0.05, transmittance);
        let mut colors = [foreground; 4];
        renderer
            .add_star_layer(&samples, 2, 2, &mut colors)
            .unwrap();
        assert_abs_diff_eq!(
            colors[0].as_vector(),
            foreground.as_vector() + transmittance * flux,
            epsilon = 1e-12
        );
        assert_eq!(colors[0].transmittance, 0.0);
        assert_eq!(&colors[1..], &[foreground; 3]);
    }
}

#[test]
fn tube_subdivision_propagates_failed_child_rays() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let good = renderer(&g, &calls, false);
    let mut samples = corners(&good, 0, 0);
    samples[0].accumulated_angular_distance = 10.0;
    let failing = renderer(&g, &calls, true);
    assert!(matches!(
        failing.handle_tube(&tube(&samples)),
        Err(RaytracerError::BelowRISCO)
    ));
}

#[test]
fn tube_rejects_degenerate_and_nonfinite_areas_at_depth_limit() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    renderer.scene.max_subdivision_depth = 0;
    let original = corners(&renderer, 0, 0);
    // Coincident sky directions, zero placeholder vectors, and invalid
    // directions must never reach an unchecked magnification division.
    for direction in [
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::zeros(),
        Vector3::new(f64::NAN, 0.0, 1.0),
        Vector3::new(f64::INFINITY, 0.0, 1.0),
    ] {
        let mut samples = original;
        for sample in &mut samples {
            sample.ray_class = RayClass::Escaped(EscapeInfo {
                x: direction.x,
                y: direction.y,
                z: direction.z,
                redshift: 1.0,
            });
        }
        assert!(matches!(
            renderer.handle_tube(&tube(&samples)),
            Err(RaytracerError::InvalidStarTube(_))
        ));
    }
    let mut samples = original;
    samples[0].ray_class = original[1].ray_class;
    assert!(matches!(
        renderer.handle_tube(&tube(&samples)),
        Err(RaytracerError::InvalidStarTube(_))
    ));
    samples = original;
    for sample in &mut samples {
        sample.ray.momentum = original[0].ray.momentum;
    }
    assert!(matches!(
        renderer.handle_tube(&tube(&samples)),
        Err(RaytracerError::InvalidStarTube(_))
    ));
    samples = original;
    samples[0].accumulated_angular_distance = f64::NAN;
    assert!(matches!(
        renderer.handle_tube(&tube(&samples)),
        Err(RaytracerError::InvalidStarTube(_))
    ));
}

#[test]
fn tube_mixed_escape_capture_subdivision_stays_finite_and_bounded() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    let mut samples = corners(&renderer, 0, 0);
    samples[0].ray_class = RayClass::Captured;
    renderer.scene.max_subdivision_depth = 2;
    calls.store(0, Ordering::Relaxed);
    let actual = renderer.handle_tube(&tube(&samples)).unwrap();
    assert!(finite_tube_flux(actual).is_ok());
    // Only the child containing the captured corner needs another split.
    assert_eq!(calls.load(Ordering::Relaxed), 10);
}

#[test]
fn tube_pass_is_skipped_without_stars_even_for_invalid_tubes() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let mut renderer = renderer(&g, &calls, false);
    let samples = corners(&renderer, 0, 0);
    let invalid = [samples[0]; 4];
    let expected = [Radiance::new(2.0, 1.0, 0.5, 0.4); 4];
    calls.store(0, Ordering::Relaxed);
    for catalogue in [None, Some(StarCatalog { stars: vec![] })] {
        renderer.scene.star_catalog = catalogue;
        let mut colors = expected;
        renderer
            .add_star_layer(&invalid, 2, 2, &mut colors)
            .unwrap();
        assert_eq!(colors, expected);
    }
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn tube_bounds_reject_unrepresentable_refinement() {
    let g = EuclideanSpace::new();
    let calls = Arc::new(AtomicUsize::new(0));
    let renderer = renderer(&g, &calls, false);
    let samples = corners(&renderer, 0, 0);
    for (top, bottom) in [
        (0.0, 0.0),
        (1.0, 0.0),
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
        (1.0, 1.0 + f64::EPSILON),
    ] {
        let mut parent = tube(&samples);
        parent.screen_bounds.top = top;
        parent.screen_bounds.bottom = bottom;
        assert!(matches!(
            renderer.trace_subdivision(&parent),
            Err(RaytracerError::InvalidStarTube(_))
        ));
    }
    assert!(SampleTube::from_buffer(&[], 0, 0, 0).is_none());
    assert!(SampleTube::from_buffer(&samples[..1], 0, 0, 1).is_none());
}
