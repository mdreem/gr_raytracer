use super::*;
use crate::geometry::euclidean::EuclideanSpace;
use crate::rendering::camera::Camera;
use crate::rendering::integrator::IntegrationConfiguration;
use crate::rendering::octree::Octree;
use crate::rendering::star_catalog::StarCatalog;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TemperatureData, TextureData, TextureMap, UVCoordinates};
use crate::scene_objects::disc::Disc;
use crate::scene_objects::objects::Objects;
use crate::scene_objects::volumetric_disc::VolumetricDisc;
use approx::assert_abs_diff_eq;
use std::sync::Arc;

struct FixedTexture(CIETristimulus);

impl TextureMap for FixedTexture {
    fn color_at_uv(
        &self,
        _: &UVCoordinates,
        _: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
        Ok(self.0)
    }
}

struct ConstantTemperature;

impl TemperatureComputer for ConstantTemperature {
    fn compute_temperature(&self, _: f64) -> Result<f64, RaytracerError> {
        Ok(1000.0)
    }
}

fn renderer<'a>(
    geometry: &'a EuclideanSpace,
    objects: Objects<'a, EuclideanSpace>,
    position: Point,
    sky: CIETristimulus,
) -> Raytracer<'a, EuclideanSpace> {
    let camera = Camera::new(
        position,
        FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
        0.7,
        3,
        3,
        0.0,
        0.0,
        0.0,
        geometry,
    )
    .unwrap();
    let scene = Scene::new(
        IntegrationConfiguration::new(100, 20.0, 0.01, 1e-8),
        objects,
        TextureData {
            celestial_map: Arc::new(FixedTexture(sky)),
        },
        geometry,
        camera,
        false,
        0.0,
    );
    Raytracer::new(scene, ToneMappingMethod::Reinhard, 1.0)
}

#[test]
fn scene_surface_and_sky_apply_straight_alpha_once() {
    let geometry = EuclideanSpace::new();
    let position = Point::new_cartesian(0.0, 2.0, 0.0, 2.0);
    let ray = Ray::new(
        1,
        1,
        position,
        FourVector::new_cartesian(-1.0, 0.0, 0.0, -1.0),
    );
    let source = Vector3::new(2.0, 1.0, 0.5);
    for opacity in [0.0, 0.1, 0.6, 1.0] {
        for sky_opacity in [0.0, 0.3, 1.0] {
            let mut objects = Objects::new(&geometry);
            objects.add_object(Box::new(Disc::new(
                0.5,
                4.0,
                Arc::new(FixedTexture(CIETristimulus::new(
                    source.x, source.y, source.z, opacity,
                ))),
                Box::new(ConstantTemperature),
                1.0,
            )));
            let sky = CIETristimulus::new(4.0, 6.0, 8.0, sky_opacity);
            let renderer = renderer(&geometry, objects, position, sky);
            let actual = renderer.scene.color_of_ray(&ray).unwrap().color;
            let expected = opacity * source + (1.0 - opacity) * sky_opacity * sky.as_vector();
            assert_abs_diff_eq!(actual.as_vector(), expected, epsilon = 1e-12);
            assert_abs_diff_eq!(
                actual.transmittance,
                (1.0 - opacity) * (1.0 - sky_opacity),
                epsilon = 1e-14
            );
        }
    }
}

#[test]
fn scene_volume_preserves_marched_light_with_or_without_a_background() {
    let geometry = EuclideanSpace::new();
    let position = Point::new_cartesian(0.0, 2.0, 0.0, 2.0);
    let ray = Ray::new(
        1,
        1,
        position,
        FourVector::new_cartesian(-1.0, 0.0, 0.0, -1.0),
    );
    let mut objects = Objects::new(&geometry);
    objects.add_object(Box::new(VolumetricDisc::new(
        1.0,
        3.0,
        Arc::new(FixedTexture(CIETristimulus::new(1.0, 1.0, 1.0, 1.0))),
        Box::new(ConstantTemperature),
        Vector3::new(0.0, 0.0, 1.0),
        4,
        42,
        500,
        0.01,
        0.5,
        10.0,
        1000.0,
        0.4,
        0.0,
        Vector3::repeat(1.0),
        1.0,
    )));
    let mut renderer = renderer(
        &geometry,
        objects,
        position,
        CIETristimulus::new(4.0, 4.0, 4.0, 1.0),
    );
    // Defer the sky exactly as catalogue rendering does, retaining gas T.
    renderer.scene.star_catalog = Some(StarCatalog {
        stars: Octree::new(vec![]),
    });
    let foreground = renderer.scene.color_of_ray(&ray).unwrap().color;
    assert!(foreground.opacity() > 0.01 && foreground.opacity() < 0.9);
    assert_abs_diff_eq!(
        foreground.as_vector(),
        Vector3::repeat(foreground.opacity()),
        epsilon = 1e-12
    );
    renderer.scene.star_catalog = None;
    let with_sky = renderer.scene.color_of_ray(&ray).unwrap().color;
    assert_abs_diff_eq!(
        with_sky.as_vector(),
        foreground.as_vector() + Vector3::repeat(4.0 * foreground.transmittance),
        epsilon = 1e-12
    );
    assert_eq!(with_sky.transmittance, 0.0);
}

struct HalfCoveredTexture;

impl TextureMap for HalfCoveredTexture {
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        _: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
        Ok(if uv.u >= 0.5 {
            CIETristimulus::new(2.0, 1.0, 0.5, 1.0)
        } else {
            CIETristimulus::new(100.0, 100.0, 100.0, 0.0)
        })
    }
}

#[test]
fn actual_supersampling_averages_covered_light_and_transmitted_sky() {
    let geometry = EuclideanSpace::new();
    let mut objects = Objects::new(&geometry);
    objects.add_object(Box::new(Disc::new(
        0.0,
        100.0,
        Arc::new(HalfCoveredTexture),
        Box::new(ConstantTemperature),
        1.0,
    )));
    let renderer = renderer(
        &geometry,
        objects,
        Point::new_cartesian(0.0, 0.0, 0.0, 10.0),
        CIETristimulus::new(20.0, 20.0, 20.0, 0.0),
    );
    let mut pixels = vec![PixelToSample {
        row: 1,
        col: 1,
        result: None,
    }];
    // Two strata on either side of the surface's x=0 coverage edge.
    renderer.supersample(2, &mut pixels).unwrap();
    let actual = pixels[0].result.unwrap();
    assert_abs_diff_eq!(
        actual.as_vector(),
        Vector3::new(1.0, 0.5, 0.25),
        epsilon = 1e-12
    );
    assert_abs_diff_eq!(actual.transmittance, 0.5, epsilon = 1e-14);
    assert_abs_diff_eq!(
        actual.over(Radiance::new(4.0, 6.0, 8.0, 0.0)).as_vector(),
        Vector3::new(3.0, 3.5, 4.25),
        epsilon = 1e-12
    );
}
