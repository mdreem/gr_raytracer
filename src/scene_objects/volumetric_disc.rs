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
use crate::scene_objects::volumetric_disc::CylinderIntersection::{
    NoIntersection, OneIntersection, Parallel, TwoIntersections,
};
use log::{error, trace};
use nalgebra::Vector3;
use noise::{NoiseFn, Perlin};

const MIN_INTERSECTION_T: f64 = 1e-9;

pub struct VolumetricDisc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
    axis: Vector3<f64>,
    e1: Vector3<f64>,
    e2: Vector3<f64>,
    perlin: Perlin,
    num_octaves: usize,
    max_steps: usize,
    step_size: f64,
    thickness: f64,
    density_multiplier: f64,
    brightness_reference_temperature: f64,
    absorption: f64,
    scattering: f64,
    noise_scale: Vector3<f64>,
    noise_offset: f64,
}

impl VolumetricDisc {
    pub fn new(
        center_disk_inner_radius: f64,
        center_disk_outer_radius: f64,
        texture_mapper: TextureMapHandle,
        temperature_computer: Box<dyn TemperatureComputer>,
        axis: Vector3<f64>,
        num_octaves: usize,
        perlin_seed: u32,
        max_steps: usize,
        step_size: f64,
        thickness: f64,
        density_multiplier: f64,
        brightness_reference_temperature: f64,
        absorption: f64,
        scattering: f64,
        noise_scale: Vector3<f64>,
        noise_offset: f64,
    ) -> Self {
        let axis = if axis.norm_squared() <= f64::EPSILON {
            Vector3::new(0.0, 0.0, 1.0)
        } else {
            axis.normalize()
        };
        let e1 = if axis.x.abs() > 0.9 {
            Vector3::new(0.0, 1.0, 0.0)
        } else {
            Vector3::new(1.0, 0.0, 0.0)
        }
        .cross(&axis)
        .normalize();
        let e2 = axis.cross(&e1).normalize();

        Self {
            center_disk_inner_radius,
            center_disk_outer_radius,
            texture_mapper,
            temperature_computer,
            axis,
            e1,
            e2,
            perlin: Perlin::new(perlin_seed),
            num_octaves,
            max_steps,
            step_size,
            thickness,
            density_multiplier,
            brightness_reference_temperature,
            absorption,
            scattering,
            noise_scale,
            noise_offset,
        }
    }

    fn compute_density(&self, p: &Vector3<f64>) -> f64 {
        let h = p.dot(&self.axis).abs();
        let r = p.cross(&self.axis).norm();

        if r <= self.center_disk_inner_radius || r >= self.center_disk_outer_radius {
            return 0.0;
        }

        // Smooth vertical falloff (Gaussian) using thickness as sigma
        let vertical_falloff = (-(h / self.thickness).powi(2)).exp();
        if vertical_falloff < 0.001 {
            return 0.0;
        }

        // Radial base density (inverse power law)
        let radial_base = (self.center_disk_inner_radius / r).powf(1.5);

        // Smooth radial boundaries
        let mut boundary_falloff = 1.0;
        boundary_falloff *= (-1.0 / (self.center_disk_outer_radius - r).powi(2).max(0.0001)).exp();
        boundary_falloff *= (-1.0 / (r - self.center_disk_inner_radius).powi(2).max(0.0001)).exp();

        // Periodic Cylindrical Noise Mapping (removes the phi-seam)
        let x_local = p.dot(&self.e1);
        let y_local = p.dot(&self.e2);
        let phi = y_local.atan2(x_local);

        // Map phi to a circle in noise space to ensure continuity
        let noise_phi_x = phi.cos() * self.noise_scale.y;
        let noise_phi_y = phi.sin() * self.noise_scale.y;

        // Use a 3D noise sample where phi-components are coordinates
        let noise_p = Vector3::new(r * self.noise_scale.x, noise_phi_x, noise_phi_y);
        let mut n = self.fbm(noise_p, 0.5);

        // Add vertical variation separately
        n += self.noise(Vector3::new(r * 0.5, h * self.noise_scale.z, phi.cos())) * 0.5;

        let n = (n + self.noise_offset).max(0.0) * self.density_multiplier;

        n * radial_base * vertical_falloff * boundary_falloff
    }

    fn get_uv(&self, p: &Vector3<f64>) -> UVCoordinates {
        let x = p.dot(&self.e1);
        let y = p.dot(&self.e2);
        let rr = (x * x + y * y).sqrt();

        let phi = y.atan2(x);
        let r = (rr - self.center_disk_inner_radius)
            / (self.center_disk_outer_radius - self.center_disk_inner_radius);

        let u = 0.5 + 0.5 * r * phi.cos();
        let v = 0.5 + 0.5 * r * phi.sin();
        UVCoordinates { u, v }
    }

    fn does_exit(&self, p: &Vector3<f64>, rd: &Vector3<f64>, t: f64) -> bool {
        let point_from = Point::new_cartesian(0.0, p[0], p[1], p[2]);
        let point_to =
            Point::new_cartesian(0.0, p[0] + rd[0] * t, p[1] + rd[1] * t, p[2] + rd[2] * t);
        if let Some(intersection) = self.intersects(&point_from, &point_to) {
            let exit = intersection.t > MIN_INTERSECTION_T;
            if exit {
                trace!(
                    "  does_exit found intersection at t={} ( > {}), breaking.",
                    intersection.t, MIN_INTERSECTION_T
                );
            }
            exit
        } else {
            false
        }
    }

    fn precompute_exit_distance(&self, ro: &Vector3<f64>, rd: &Vector3<f64>) -> Option<f64> {
        let max_distance = self.step_size * self.max_steps as f64;
        let to = ro + rd * max_distance;
        let intersection = self.intersects_cylinder(
            ro,
            &to,
            self.center_disk_inner_radius,
            self.center_disk_outer_radius,
            &self.axis,
        );

        match intersection {
            NoIntersection | Parallel => None,
            OneIntersection(t) => (t > MIN_INTERSECTION_T).then_some(t * max_distance),
            TwoIntersections(t1, t2) => {
                if t1 > MIN_INTERSECTION_T {
                    Some(t1 * max_distance)
                } else if t2 > MIN_INTERSECTION_T {
                    Some(t2 * max_distance)
                } else {
                    None
                }
            }
        }
    }

    // https://www.scratchapixel.com/lessons/3d-basic-rendering/volume-rendering-for-developers/ray-marching-get-it-right.html
    fn raymarch_constant_step(
        &self,
        ro: &Vector3<f64>,
        rd: &Vector3<f64>,
        redshift: f64,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.raymarch_constant_step_internal(ro, rd, redshift, true)
    }

    fn raymarch_constant_step_internal(
        &self,
        ro: &Vector3<f64>,
        rd: &Vector3<f64>,
        redshift: f64,
        use_cached_exit: bool,
    ) -> Result<CIETristimulus, RaytracerError> {
        let sigma_a = self.absorption;
        let sigma_s = self.scattering;

        let mut accum_color = CIETristimulus::new(0.0, 0.0, 0.0, 0.0);
        let mut transparency = 1.0;
        let mut alpha_weighted_sum = 0.0;
        let mut alpha_weight_total = 0.0;

        let d_s = self.step_size;
        let mut step_count = 0;
        let mut d_o = 0.0;
        let cached_exit_distance = if use_cached_exit {
            self.precompute_exit_distance(ro, rd)
        } else {
            None
        };

        for i in 0..self.max_steps {
            let p = ro + rd * d_o;
            d_o += d_s;

            step_count += 1;
            let density = self.compute_density(&p);

            if density > 0.0 {
                let r_dist = p.cross(&self.axis).norm();
                let temperature = self.temperature_computer.compute_temperature(r_dist)?;
                let uv = self.get_uv(&p);
                let light_color = self.texture_mapper.color_at_uv(
                    &uv,
                    &crate::rendering::texture::TemperatureData {
                        temperature,
                        redshift,
                    },
                );

                let sample_attenuation = (-d_s * density * (sigma_a + sigma_s)).exp();
                transparency *= sample_attenuation;

                let travel_density = d_s; // Ignore the light ray here for now.
                let light_attenuation = (-density * travel_density * (sigma_a + sigma_s)).exp();

                // Stefan-Boltzmann law: emission intensity scales with T^4.
                // Use a reference temperature for normalization to boost brightness.
                let intensity_factor =
                    (temperature / self.brightness_reference_temperature).powi(4);

                let emission_weight = transparency * light_attenuation * sigma_s * density * d_s;
                let step_emission = light_color.mul_color_part(emission_weight * intensity_factor);

                // Keep texture alpha influence separate from per-step emission accumulation.
                let alpha_sample_weight = density * d_s;
                alpha_weighted_sum += light_color.alpha.clamp(0.0, 1.0) * alpha_sample_weight;
                alpha_weight_total += alpha_sample_weight;

                accum_color.x += step_emission.x;
                accum_color.y += step_emission.y;
                accum_color.z += step_emission.z;
            }

            // Use precomputed exit distance when available (fast path), fallback to legacy check.
            let exited = if let Some(exit_distance) = cached_exit_distance {
                d_o >= exit_distance
            } else {
                self.does_exit(&p, &rd, d_s)
            };
            if exited {
                trace!("  Ray exited disc at step {}, position {:?}", i, p,);
                break;
            }
        }

        // Final alpha combines physical opacity with texture alpha, applied once at the end.
        let physical_opacity = 1.0 - transparency;
        let texture_alpha = if alpha_weight_total > 0.0 {
            alpha_weighted_sum / alpha_weight_total
        } else {
            1.0
        };
        accum_color.alpha = physical_opacity * texture_alpha;

        trace!(
            "  Computed from {:?} to {:?} with step_count={}",
            ro,
            ro + rd * d_o,
            step_count
        );
        trace!("  resulting color: {:?}", accum_color);
        Ok(accum_color)
    }

    fn fbm(&self, x: Vector3<f64>, h: f64) -> f64 {
        let g = (-h).exp2();
        let mut frequency = 4.0;
        let mut amplitude = 1.0;
        let mut t = 0.0;

        for _ in 0..self.num_octaves {
            t += amplitude * self.noise(x * frequency);
            frequency *= 2.0;
            amplitude *= g;
        }
        t
    }

    fn noise(&self, p: Vector3<f64>) -> f64 {
        self.perlin.get([p[0], p[1], p[2]])
    }

    fn intersects_clipped_cylinder(
        &self,
        from: &Vector3<f64>,
        to: &Vector3<f64>,
        radius: f64,
        axis: &Vector3<f64>,
        half_height: f64,
    ) -> CylinderIntersection {
        let segment_vector = to - from;
        let segment_length = segment_vector.norm();

        if segment_length < 1e-12 {
            return NoIntersection;
        }
        let d = segment_vector / segment_length;

        let v = from.cross(axis);
        let w = d.cross(axis);

        let a = w.dot(&w);
        let b = 2.0 * v.dot(&w);
        let c = v.dot(&v) - radius * radius;

        if a < 1e-10 {
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
        &self,
        from: &Vector3<f64>,
        to: &Vector3<f64>,
        radius: f64,
        axis: &Vector3<f64>,
        pos: f64,
    ) -> Option<f64> {
        let segment_vector = to - from;
        let segment_length = segment_vector.norm();
        if segment_length < 1e-12 {
            return None;
        }
        let n = segment_vector / segment_length;

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
        &self,
        from: &Vector3<f64>,
        to: &Vector3<f64>,
        inner_radius: f64,
        outer_radius: f64,
        axis: &Vector3<f64>,
    ) -> CylinderIntersection {
        let mut hits = Vec::new();
        let dir = to - from;

        // Use 3.0 * thickness to capture the Gaussian tail
        let capture_height = self.thickness * 3.0;

        // Check Outer Tube
        match self.intersects_clipped_cylinder(from, to, outer_radius, axis, capture_height) {
            OneIntersection(t) => hits.push(t),
            TwoIntersections(t1, t2) => {
                hits.push(t1);
                hits.push(t2);
            }
            _ => {}
        }
        // Check Inner Tube
        match self.intersects_clipped_cylinder(from, to, inner_radius, axis, capture_height) {
            OneIntersection(t) => hits.push(t),
            TwoIntersections(t1, t2) => {
                hits.push(t1);
                hits.push(t2);
            }
            _ => {}
        }
        // Check Caps
        for pos in [capture_height, -capture_height] {
            if let Some(t) = self.intersects_cap(from, to, outer_radius, axis, pos) {
                let p = from + t * dir;
                let r_sq = p.cross(axis).norm_squared();
                if r_sq >= inner_radius * inner_radius {
                    hits.push(t);
                }
            }
        }

        hits.sort_by(|a, b| a.total_cmp(b));

        if hits.is_empty() {
            NoIntersection
        } else if hits.len() == 1 {
            OneIntersection(hits[0])
        } else {
            TwoIntersections(hits[0], hits[1])
        }
    }
}

#[derive(Debug)]
pub enum CylinderIntersection {
    NoIntersection,
    Parallel,
    OneIntersection(f64),
    TwoIntersections(f64, f64),
}

impl Hittable for VolumetricDisc {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection> {
        let y_start_spatial = y_start.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.get_spatial_vector_cartesian();

        let direction = y_end_spatial - y_start_spatial;

        let cylinder_intersection = self.intersects_cylinder(
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
            NoIntersection => return None,
            Parallel => return None,
            OneIntersection(t) => {
                if t > MIN_INTERSECTION_T {
                    t
                } else {
                    return None;
                }
            }
            TwoIntersections(t1, t2) => {
                if t1 > MIN_INTERSECTION_T {
                    t1
                } else if t2 > MIN_INTERSECTION_T {
                    t2
                } else {
                    return None;
                }
            }
        };

        let intersection_point = y_start_spatial + t * direction;

        if !(0.0..=1.0).contains(&t) {
            error!("VolumetricDisc: Intersection t={} out of bounds.", t);
            return None;
        }

        trace!(
            "VolumetricDisc intersection at t = {}, point = {:?} with distance {}",
            t,
            intersection_point,
            intersection_point.norm()
        );

        let uv = self.get_uv(&intersection_point);

        Some(Intersection {
            uv,
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

        self.raymarch_constant_step(&ro, &rd, color_computation_data.temperature_data.redshift)
            .unwrap_or_else(|e| {
                log::warn!("Raymarching failed: {}", e);
                CIETristimulus::new(0.0, 0.0, 0.0, 0.0)
            })
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

    fn temperature_of_emitter(
        &self,
        point: &Point,
        geometry: &dyn Geometry,
    ) -> Result<f64, RaytracerError> {
        let r = geometry.get_radial_coordinate(point);
        let temperature = self.temperature_computer.compute_temperature(r)?;
        Ok(temperature)
    }
}

impl SceneObject for VolumetricDisc {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rendering::color::CIETristimulusNormalization::NoNormalization;
    use crate::rendering::color::Color;
    use crate::rendering::texture::{CheckerMapper, TemperatureData, TextureMap};
    use approx::assert_abs_diff_eq;
    use std::sync::Arc;

    struct DummyTemperatureComputer;

    impl TemperatureComputer for DummyTemperatureComputer {
        fn compute_temperature(&self, _r: f64) -> Result<f64, RaytracerError> {
            Ok(1000.0)
        }
    }

    struct FixedTextureMap {
        color: CIETristimulus,
    }

    impl TextureMap for FixedTextureMap {
        fn color_at_uv(
            &self,
            _uv: &UVCoordinates,
            _temperature_data: &TemperatureData,
        ) -> CIETristimulus {
            self.color
        }
    }

    fn create_disc() -> VolumetricDisc {
        VolumetricDisc::new(
            1.0,
            3.0,
            Arc::new(CheckerMapper::new(
                3.0,
                5.0,
                5.0,
                Color::new(255, 0, 0, 255),
                Color::new(0, 0, 255, 255),
                NoNormalization,
            )),
            Box::new(DummyTemperatureComputer),
            Vector3::new(0.0, 0.0, 1.0),
            4,
            42,
            500,
            0.01,
            0.5,
            10.0,
            1000.0,
            0.2,
            0.2,
            Vector3::new(1.0, 1.0, 1.0),
            1.0,
        )
    }

    #[test]
    fn test_volumetric_disc_intersection_exists() {
        let disc = create_disc();
        let y_start = Point::new_cartesian(0.0, 0.5, 0.0, 0.0);
        let y_end = Point::new_cartesian(0.0, 1.5, 0.0, 0.0);

        assert!(disc.intersects(&y_start, &y_end).is_some());
    }

    #[test]
    fn test_volumetric_disc_intersection_miss_above_caps() {
        let disc = create_disc();
        let y_start = Point::new_cartesian(0.0, 1.5, 0.0, 2.0);
        let y_end = Point::new_cartesian(0.0, 2.5, 0.0, 2.0);

        assert!(disc.intersects(&y_start, &y_end).is_none());
    }

    #[test]
    fn test_volumetric_disc_density_inside_and_outside() {
        let disc = create_disc();

        let inside = Vector3::new(2.0, 0.0, 0.0);
        let outside_radial = Vector3::new(0.2, 0.0, 0.0);
        let outside_vertical = Vector3::new(2.0, 0.0, 2.0);

        assert!(disc.compute_density(&inside) > 0.0);
        assert_eq!(disc.compute_density(&outside_radial), 0.0);
        assert_eq!(disc.compute_density(&outside_vertical), 0.0);
    }

    #[test]
    fn test_volumetric_disc_raymarch_produces_opacity() {
        let disc = create_disc();
        let ro = Vector3::new(2.0, 0.0, 0.0);
        let rd = Vector3::new(1.0, 0.0, 0.0);

        let color = disc
            .raymarch_constant_step(&ro, &rd, 1.0)
            .expect("raymarch should succeed");

        assert!(color.alpha > 0.0);
        assert!(color.alpha <= 1.0);
    }

    #[test]
    fn test_volumetric_disc_raymarch_cached_exit_matches_legacy() {
        let disc = VolumetricDisc::new(
            1.0,
            3.0,
            Arc::new(FixedTextureMap {
                color: CIETristimulus::new(2.0, 1.0, 0.5, 0.7),
            }),
            Box::new(DummyTemperatureComputer),
            Vector3::new(0.0, 0.0, 1.0),
            4,
            42,
            500,
            0.01,
            0.5,
            10.0,
            1000.0,
            0.2,
            0.2,
            Vector3::new(1.0, 1.0, 1.0),
            1.0,
        );
        let ro = Vector3::new(2.0, 0.0, 0.0);
        let rd = Vector3::new(1.0, 0.0, 0.0);

        let cached = disc
            .raymarch_constant_step_internal(&ro, &rd, 1.0, true)
            .expect("raymarch should succeed");
        let legacy = disc
            .raymarch_constant_step_internal(&ro, &rd, 1.0, false)
            .expect("raymarch should succeed");

        assert_abs_diff_eq!(cached.x, legacy.x, epsilon = 1e-10);
        assert_abs_diff_eq!(cached.y, legacy.y, epsilon = 1e-10);
        assert_abs_diff_eq!(cached.z, legacy.z, epsilon = 1e-10);
        assert_abs_diff_eq!(cached.alpha, legacy.alpha, epsilon = 1e-10);
    }
}
