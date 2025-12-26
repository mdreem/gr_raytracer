use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use log::info;
use nalgebra::{Vector3, inf};
use noise::{NoiseFn, Perlin};

const NUM_OCTAVES: usize = 8;
const MAX_STEPS: usize = 1000;

pub struct VolumetricDisc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
}

impl VolumetricDisc {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: TextureMapHandle,
        temperature_computer: Box<dyn TemperatureComputer>,
    ) -> Self {
        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
            temperature_computer,
        }
    }

    fn raymarch_constant_step(&self, ro: &Vector3<f64>, rd: &Vector3<f64>) -> CIETristimulus {
        let sigma_a = 0.5; // absorption coefficient
        let sigma_s = 0.8; // scattering coefficient

        let background_color = CIETristimulus::new(0.0, 0.0, 0.0, 0.0);
        let light_color = CIETristimulus::new(1.0, 1.0, 1.0, 1.0);
        let mut result = CIETristimulus::new(0.0, 0.0, 0.0, 0.0);

        let mut transparency = 1.0;

        let max_distance = 10.0 * self.center_disk_outer_radius;
        let d_s = max_distance / (MAX_STEPS as f64);

        // TODO count number of times we went through the middle.

        let mut hit_ring = false;
        let mut exited_ring = false;

        let mut step_count = 0;

        let mut dO = 0.0;
        for i in 0..MAX_STEPS {
            let p = ro + rd * dO;
            dO += d_s;

            step_count += 1;

            let length = p.norm();
            // TODO: ensure intersection is always inside the ring first.
            if length < self.center_disk_inner_radius - 1e-2
                || length > self.center_disk_outer_radius
            {
                if !exited_ring {
                    info!(
                        "  Exited ring at step {}, position {:?}, length {}",
                        i, p, length
                    );
                    // Only consider entering along the direction of the extended side.
                    if i == 0 {
                        return CIETristimulus::new(0.0, 0.0, 0.0, 0.0);
                    }
                }

                exited_ring = true;
                break;
            }
            if exited_ring {
                info!(
                    "  Re-entered ring at step {}, position {:?}, length {}",
                    i, p, length
                );
                exited_ring = false;
            }

            hit_ring = true;

            let n = fbm(p, 0.5);
            let n = n * 0.5 + 0.5;

            let mut density = n;
            let disc_thickness = (-p[2].abs() / 0.1).exp(); // TODO: make the denominator a parameter.
            density *= disc_thickness;

            let sample_attenuation = (-d_s * density * (sigma_a + sigma_s)).exp();
            transparency *= sample_attenuation;

            if density > 0.0 {
                let light_attenuation = (-density * (sigma_a + sigma_s)).exp();
                result = result
                    + ((transparency * light_attenuation * sigma_s * density * d_s) * light_color);
            }
        }
        info!(
            "  Computed from {:?} to {:?} with step_count={}",
            ro,
            ro + rd * dO,
            step_count
        );
        let mut res = (transparency * background_color) + result;
        info!("  resulting color: {:?}", res);
        res
    }
}

/// Wrapper for Perlin noise function.
fn noise(p: Vector3<f64>) -> f64 {
    let perlin = Perlin::new(1);
    let val = perlin.get([p[0], p[1], p[2]]);
    val
}

/// Fractional Brownian Motion (fbm) using `Vector3<f64>` positions.
/// `h` is the Hölder exponent controlling amplitude falloff.
fn fbm(x: Vector3<f64>, h: f64) -> f64 {
    // gain per octave: 2^{-h}
    let g = (-h).exp2();
    let mut frequency = 1.0_f64;
    let mut amplitude = 1.0_f64;
    let mut t = 0.0_f64;

    for _ in 0..NUM_OCTAVES {
        t += amplitude * noise(x * frequency);
        frequency *= 2.0;
        amplitude *= g;
    }
    t
}

impl Hittable for VolumetricDisc {
    // TODO: ensure sphere is not skipped due to large step sizes.
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        let y_start_spatial = y_start.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        let n_start = y_start_spatial.norm_squared().sqrt();
        let n_end = y_end_spatial.norm_squared().sqrt();

        // If both points are inside inner radius, no intersection.
        if n_end < self.center_disk_inner_radius && n_start < self.center_disk_inner_radius {
            return None;
        }
        // If both points are outside outer radius, no intersection.
        if n_start > self.center_disk_outer_radius && n_end > self.center_disk_outer_radius {
            return None;
        }

        // Check if the segment crosses the outer radius.
        let t = if (n_start < self.center_disk_outer_radius
            && n_end > self.center_disk_outer_radius)
            || (n_start > self.center_disk_outer_radius && n_end < self.center_disk_outer_radius)
        {
            // crossed outer radius
            (self.center_disk_outer_radius - n_start) / (n_end - n_start)
        }
        // Both inside outer radius, check inner radius.
        else if (n_start < self.center_disk_inner_radius && n_end > self.center_disk_inner_radius)
            || (n_start > self.center_disk_inner_radius && n_end < self.center_disk_inner_radius)
        {
            (self.center_disk_inner_radius - n_start) / (n_end - n_start)
        } else {
            return None;
        };

        let intersection_point = y_start_spatial + t * direction;

        info!(
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
