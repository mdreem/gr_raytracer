use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::{CIETristimulus, Color, srgb_to_xyz};
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use crate::scene_objects::volumetric_disc::CylinderIntersection::{
    NoIntersection, OneIntersection, Parallel, TwoIntersections,
};
use log::{info, trace};
use nalgebra::{DimMul, Vector3, inf};
use noise::{NoiseFn, Perlin};

const NUM_OCTAVES: usize = 8;
const MAX_STEPS: usize = 10000;
const STEP_SIZE: f64 = 0.01;
const DISC_THICKNESS: f64 = 0.02;

pub struct VolumetricDisc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
    axis: Vector3<f64>,
    perlin: Perlin,
}

impl VolumetricDisc {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: TextureMapHandle,
        temperature_computer: Box<dyn TemperatureComputer>,
        axis: Vector3<f64>,
    ) -> Self {
        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
            temperature_computer,
            axis: axis.normalize(),
            perlin: Perlin::new(1),
        }
    }

    fn compute_density(&self, p: &Vector3<f64>) -> f64 {
        let h = p.dot(&self.axis).abs();
        let r = p.cross(&self.axis).norm();

        if r <= self.center_disk_inner_radius
            || r >= self.center_disk_outer_radius
            || h > DISC_THICKNESS
        {
            return 0.0;
        }

        let mut thickness = 1.0;

        thickness *= (-1.0 / (self.center_disk_outer_radius - r).powi(4).max(0.00001)).exp();
        thickness *= (-1.0 / (r - self.center_disk_inner_radius).powi(4).max(0.00001)).exp();

        let n = self.fbm(*p, 0.5);
        let n = 10.55 * (n + 0.90);
        n * thickness
    }

    fn does_exit(&self, p: &Vector3<f64>, rd: &Vector3<f64>, t: f64) -> bool {
        let point_from = Point::new_cartesian(0.0, p[0], p[1], p[2]);
        let point_to =
            Point::new_cartesian(0.0, p[0] + rd[0] * t, p[1] + rd[1] * t, p[2] + rd[2] * t);
        if let Some(intersection) = self.intersects(&point_from, &point_to) {
            let exit = intersection.t > 1e-9;
            if exit {
                trace!(
                    "  does_exit found intersection at t={} ( > 1e-9), breaking.",
                    intersection.t
                );
            }
            exit
        } else {
            false
        }
    }

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/volume-rendering-for-developers/ray-marching-get-it-right.html
    fn raymarch_constant_step(&self, ro: &Vector3<f64>, rd: &Vector3<f64>) -> CIETristimulus {
        let sigma_a = 0.9; // absorption coefficient
        let sigma_s = 0.9; // scattering coefficient

        let light_color = srgb_to_xyz(&Color::new(255, 100, 0, 255)); // orange light

        let mut result = 0.0;
        let mut transparency = 1.0;

        let d_s = STEP_SIZE;
        let mut step_count = 0;
        let mut dO = 0.0;

        for i in 0..MAX_STEPS {
            let p = ro + rd * dO;
            dO += d_s;

            step_count += 1;
            let density = self.compute_density(&p);

            let sample_attenuation = (-d_s * density * (sigma_a + sigma_s)).exp();
            transparency *= sample_attenuation;

            let travel_density = d_s; // Ignore the light ray here for now.
            let light_attenuation = (-density * travel_density * (sigma_a + sigma_s)).exp();
            result += transparency * light_attenuation * sigma_s * density * d_s;

            // Check if ray will exist the disc in the next step.
            if self.does_exit(&p, &rd, d_s) {
                trace!("  Ray exited disc at step {}, position {:?}", i, p,);
                break;
            }
        }
        trace!(
            "  Computed from {:?} to {:?} with step_count={}",
            ro,
            ro + rd * dO,
            step_count
        );
        let res = light_color.mul_all_parts(result);
        trace!("  resulting color: {:?}", res);
        res
    }

    fn fbm(&self, x: Vector3<f64>, h: f64) -> f64 {
        let g = (-h).exp2();
        let mut frequency = 4.0;
        let mut amplitude = 1.0;
        let mut t = 0.0;

        for _ in 0..NUM_OCTAVES {
            t += amplitude * self.noise(x * frequency);
            frequency *= 2.0;
            amplitude *= g;
        }
        t
    }

    fn noise(&self, p: Vector3<f64>) -> f64 {
        self.perlin.get([p[0], p[1], p[2]])
    }
}

#[derive(Debug)]
enum CylinderIntersection {
    NoIntersection,
    Parallel,
    OneIntersection(f64),
    TwoIntersections(f64, f64),
}

fn intersects_clipped_cylinder(
    from: &Vector3<f64>,
    to: &Vector3<f64>,
    radius: f64,
    axis: &Vector3<f64>,
    half_height: f64,
) -> CylinderIntersection {
    let d = (to - from).normalize();
    let segment_length = (to - from).norm();

    if segment_length < 1e-12 {
        return NoIntersection;
    }

    // Ray: P(t) = F + tD. Distance to axis A through origin: ||(F + tD) x A||^2 = R^2
    // Let V = F x A, W = D x A
    // ||V||^2 + 2t(V.W) + t^2||W||^2 = R^2
    let v = from.cross(axis);
    let w = d.cross(axis);

    let a = w.dot(&w);
    let b = 2.0 * v.dot(&w);
    let c = v.dot(&v) - radius * radius;

    if a < 1e-10 {
        // Ray parallel to axis
        if v.norm_squared() > radius * radius {
            return NoIntersection;
        }
        return Parallel;
    }

    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return NoIntersection;
    }

    let sqrt_disc = discriminant.sqrt();
    let t1_dist = (-b - sqrt_disc) / (2.0 * a);
    let t2_dist = (-b + sqrt_disc) / (2.0 * a);

    let mut hits = Vec::new();
    for dist in [t1_dist, t2_dist] {
        let t = dist / segment_length;
        if (0.0..=1.0).contains(&t) {
            let p = from + t * (to - from);
            if p.dot(axis).abs() <= half_height {
                hits.push(t);
            }
        }
    }

    if hits.is_empty() {
        return NoIntersection;
    }
    if hits.len() == 1 {
        return OneIntersection(hits[0]);
    }
    TwoIntersections(hits[0].min(hits[1]), hits[0].max(hits[1]))
}

fn intersects_cap(
    from: &Vector3<f64>,
    to: &Vector3<f64>,
    radius: f64,
    axis: &Vector3<f64>,
    pos: f64,
) -> Option<f64> {
    let n = (to - from).normalize();
    let segment_vector = to - from;

    let denom = n.dot(axis);
    if denom.abs() < 1e-10 {
        return None;
    }

    let t = (pos - from.dot(axis)) / segment_vector.dot(axis);

    if !(0.0..=1.0).contains(&t) {
        return None;
    }

    let p = from + t * segment_vector;
    let radial_dist_sq = p.cross(axis).norm_squared();
    if radial_dist_sq > radius * radius {
        return None;
    }

    Some(t)
}

fn intersects_cylinder(
    from: &Vector3<f64>,
    to: &Vector3<f64>,
    inner_radius: f64,
    outer_radius: f64,
    axis: &Vector3<f64>,
) -> CylinderIntersection {
    let mut hits = Vec::new();
    let dir = to - from;

    // 1. Check Outer Tube
    match intersects_clipped_cylinder(from, to, outer_radius, axis, DISC_THICKNESS) {
        OneIntersection(t) => hits.push(t),
        TwoIntersections(t1, t2) => {
            hits.push(t1);
            hits.push(t2);
        }
        _ => {}
    }
    // 2. Check Inner Tube
    match intersects_clipped_cylinder(from, to, inner_radius, axis, DISC_THICKNESS) {
        OneIntersection(t) => hits.push(t),
        TwoIntersections(t1, t2) => {
            hits.push(t1);
            hits.push(t2);
        }
        _ => {}
    }
    // 3. Check Caps
    for pos in [DISC_THICKNESS, -DISC_THICKNESS] {
        if let Some(t) = intersects_cap(from, to, outer_radius, axis, pos) {
            let p = from + t * dir;
            let r_sq = p.cross(axis).norm_squared();
            if r_sq >= inner_radius * inner_radius {
                hits.push(t);
            }
        }
    }

    hits.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if hits.is_empty() {
        NoIntersection
    } else if hits.len() == 1 {
        OneIntersection(hits[0])
    } else {
        TwoIntersections(hits[0], hits[1])
    }
}

impl Hittable for VolumetricDisc {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        let y_start_spatial = y_start.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.get_spatial_vector_cartesian();

        let direction = y_end_spatial - y_start_spatial;

        let cylinder_intersection = intersects_cylinder(
            &y_start_spatial,
            &y_end_spatial,
            self.center_disk_inner_radius,
            self.center_disk_outer_radius,
            &self.axis,
        );

        match &cylinder_intersection {
            NoIntersection => {}
            _ => {
                trace!(
                    "VolumetricDisc: Cylinder intersection result: {:?}",
                    cylinder_intersection
                );
            }
        }

        let t = match cylinder_intersection {
            NoIntersection => {
                return None;
            }
            Parallel => {
                return None;
            }
            OneIntersection(t) => t,
            TwoIntersections(t1, _t2) => t1, // Always take the first hit
        };

        let intersection_point = y_start_spatial + t * direction;

        if !(0.0..=1.0).contains(&t) {
            panic!("VolumetricDisc: Intersection t={} out of bounds.", t);
        }

        trace!(
            "VolumetricDisc intersection at t = {}, point = {:?} with distance {}",
            t,
            intersection_point,
            intersection_point.norm()
        );

        Some(Intersection {
            uv: UVCoordinates { u: 0.0, v: 0.0 },
            intersection_point: Point::new_cartesian(
                0.0,
                intersection_point[0],
                intersection_point[1],
                intersection_point[2],
            ),
            t,
            direction: FourVector::new_cartesian(0.0, direction[0], direction[1], direction[2]),
        })
    }

    fn color_at_uv(&self, color_computation_data: &ColorComputationData) -> CIETristimulus {
        let ro = color_computation_data
            .intersection_point
            .get_spatial_vector_cartesian();
        let rd = Vector3::<f64>::new(
            color_computation_data.direction[1],
            color_computation_data.direction[2],
            color_computation_data.direction[3],
        )
        .normalize();

        self.raymarch_constant_step(&ro, &rd)
    }

    fn energy_of_emitter(
        &self,
        geometry: &dyn Geometry,
        step: &Step,
    ) -> Result<f64, RaytracerError> {
        let position = step.x;
        let velocity = geometry.get_circular_orbit_velocity_at(&position)?;
        let momentum = step.p;
        Ok(geometry.inner_product(&position, &velocity, &momentum))
    }

    fn temperature_of_emitter(&self, point: &Point) -> Result<f64, RaytracerError> {
        let r = point.get_as_spherical()[0];
        let temperature = self.temperature_computer.compute_temperature(r)?;
        Ok(temperature)
    }
}

impl SceneObject for VolumetricDisc {}

mod tests {
    use crate::rendering::temperature::TemperatureComputer;
    use crate::rendering::texture::TextureMapHandle;
    use crate::scene_objects::volumetric_disc::VolumetricDisc;
    use nalgebra::Vector3;

    struct DummyTemperatureComputer;
    impl TemperatureComputer for DummyTemperatureComputer {
        fn compute_temperature(
            &self,
            _r: f64,
        ) -> Result<f64, crate::rendering::raytracer::RaytracerError> {
            Ok(1000.0)
        }
    }

    #[test]
    fn test_fbm() {
        let disc = VolumetricDisc::new(
            4.05,
            11.0,
            TextureMapHandle::new_checker(
                3.0,
                5.0,
                5.0,
                crate::rendering::color::Color::new(255, 0, 0, 255),
                crate::rendering::color::Color::new(0, 0, 255, 255),
                crate::rendering::color::CIETristimulusNormalization::NoNormalization,
            ),
            Box::new(DummyTemperatureComputer),
            Vector3::new(0.0, 0.0, 1.0),
        );

        for x in 0..10 {
            for y in 0..10 {
                for z in 0..10 {
                    let p = nalgebra::Vector3::new(x as f64 * 0.1, y as f64 * 0.1, z as f64 * 0.1);
                    let value = disc.fbm(p, 0.5);
                    println!("fbm value at {:?} is {}", p, value);
                }
            }
        }
        let p = nalgebra::Vector3::new(1.0, 2.0, 13.0);
        let value = disc.fbm(p, 0.5);
        println!("fbm value at {:?} is {}", p, value);
        // assert!(false);
    }
}
