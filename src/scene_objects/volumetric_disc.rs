use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::redshift::RayFrequencyData;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection, Segment};
use crate::scene_objects::objects::SceneObject;
use crate::scene_objects::volumetric_disc::CylinderIntersection::{
    NoIntersection, OneIntersection, Parallel, TwoIntersections,
};
use log::{error, trace};
use nalgebra::Vector3;
use noise::{NoiseFn, Perlin};

const MIN_INTERSECTION_T: f64 = 1e-9;
const CAPTURE_HEIGHT_FACTOR: f64 = 3.0;
/// Stop marching once accumulated transparency drops below this: everything
/// farther along the ray is occluded. The residual error is bounded by
/// transparency times the radiance behind it; this is an HDR renderer, so
/// radiance is unbounded (boosted inner-disc regions reach ~1e5 linear) and
/// the threshold is chosen conservatively small rather than the classic 1e-3.
const TRANSPARENCY_EARLY_EXIT: f64 = 1e-5;
/// `-ln(0.001)` = `ln(1000)`, where 0.001 is the vertical falloff below which
/// a sample counts as outside the disc. The support test compares the
/// falloff's exponent against this rather than evaluating `exp` on every
/// sample: `exp(-q) >= f` iff `q <= -ln(f)`.
const MAXIMUM_VERTICAL_FALLOFF_EXPONENT: f64 = 6.907_755_278_982_137;

/// A point expressed in the disc's cylindrical frame.
struct CylindricalFrame {
    /// Distance from the disc's mid-plane.
    h: f64,
    /// Cylindrical radius about the disc axis.
    r: f64,
    cos_phi: f64,
    sin_phi: f64,
}

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
    step_size: f64,
    thickness: f64,
    density_multiplier: f64,
    brightness_reference_temperature: f64,
    absorption: f64,
    scattering: f64,
    noise_scale: Vector3<f64>,
    noise_offset: f64,
}

struct SegmentState {
    distance_accumulated: f64,
    transparency: f64,
    alpha_weighted_sum: f64,
    alpha_weight_total: f64,
    color_accumulated: CIETristimulus,
}

#[derive(PartialEq)]
enum RaymarchResult {
    Continue,
    AlreadyOpaque,
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
        // Retained for configuration compatibility; the march is bounded
        // by the ray's step slice since plan 02, not by a step budget.
        _max_steps: usize,
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

    /// A point resolved into the disc's own cylindrical frame: height above
    /// the mid-plane, cylindrical radius, and the cosine/sine of the azimuth.
    ///
    /// `(e1, e2, axis)` is orthonormal, so the cylindrical radius is just the
    /// length of the `(e1, e2)` projection and the azimuth's cosine and sine
    /// are that projection normalized - no cross product, no `atan2`, and no
    /// `sin`/`cos` to undo it. Every marched sample goes through here.
    fn cylindrical_frame(&self, p: &Vector3<f64>) -> CylindricalFrame {
        let x_local = p.dot(&self.e1);
        let y_local = p.dot(&self.e2);
        // `sqrt(x^2 + y^2)`, not `hypot`: glibc's `hypot` pays for extreme-range
        // exactness that these ordinary magnitudes never need, and it showed up
        // as ~2% of a volumetric render on its own.
        let r = (x_local * x_local + y_local * y_local).sqrt();
        let (cos_phi, sin_phi) = if r > 0.0 {
            (x_local / r, y_local / r)
        } else {
            (1.0, 0.0)
        };
        CylindricalFrame {
            h: p.dot(&self.axis).abs(),
            r,
            cos_phi,
            sin_phi,
        }
    }

    /// Cheap bounding test for whether a point is inside the disc's support region.
    fn is_inside_disc(&self, frame: &CylindricalFrame) -> bool {
        if frame.r <= self.center_disk_inner_radius || frame.r >= self.center_disk_outer_radius {
            return false;
        }

        // `vertical_falloff >= MINIMUM_VERTICAL_FALLOFF` written without the
        // exponential: exp(-q) >= f is q <= -ln(f). Keeps the rejection path
        // (most samples) free of a transcendental call.
        let q = frame.h / self.thickness;
        q * q <= MAXIMUM_VERTICAL_FALLOFF_EXPONENT
    }

    // Smooth vertical falloff (Gaussian) using thickness as sigma
    fn compute_vertical_falloff(&self, h: f64) -> f64 {
        let q = h / self.thickness;
        (-(q * q)).exp()
    }

    #[allow(dead_code)] // Used by this module's tests.
    fn compute_density(&self, p: &Vector3<f64>) -> f64 {
        self.compute_density_in_frame(&self.cylindrical_frame(p))
    }

    fn compute_density_in_frame(&self, frame: &CylindricalFrame) -> f64 {
        if !self.is_inside_disc(frame) {
            return 0.0;
        }

        let (h, r) = (frame.h, frame.r);

        // Part of the density profile, not only of the support test above.
        let vertical_falloff = self.compute_vertical_falloff(h);

        // Radial base density (inverse power law). x^1.5 = x * sqrt(x), which
        // is a multiply and a hardware square root instead of a `pow` call.
        let radial_ratio = self.center_disk_inner_radius / r;
        let radial_base = radial_ratio * radial_ratio.sqrt();

        // Smooth radial boundaries
        let outer_gap = self.center_disk_outer_radius - r;
        let inner_gap = r - self.center_disk_inner_radius;
        let mut boundary_falloff = 1.0;
        boundary_falloff *= (-1.0 / (outer_gap * outer_gap).max(0.0001)).exp();
        boundary_falloff *= (-1.0 / (inner_gap * inner_gap).max(0.0001)).exp();

        // Periodic Cylindrical Noise Mapping (removes the phi-seam):
        // map phi to a circle in noise space to ensure continuity.
        let noise_phi_x = frame.cos_phi * self.noise_scale.y;
        let noise_phi_y = frame.sin_phi * self.noise_scale.y;

        // Use a 3D noise sample where phi-components are coordinates
        let noise_p = Vector3::new(r * self.noise_scale.x, noise_phi_x, noise_phi_y);
        let mut n = self.fbm(noise_p, 0.5);

        // Add vertical variation separately
        n += self.noise(Vector3::new(r * 0.5, h * self.noise_scale.z, frame.cos_phi)) * 0.5;

        let n = (n + self.noise_offset).max(0.0) * self.density_multiplier;

        n * radial_base * vertical_falloff * boundary_falloff
    }

    fn get_uv(&self, p: &Vector3<f64>) -> UVCoordinates {
        self.get_uv_in_frame(&self.cylindrical_frame(p))
    }

    fn get_uv_in_frame(&self, frame: &CylindricalFrame) -> UVCoordinates {
        let r = (frame.r - self.center_disk_inner_radius)
            / (self.center_disk_outer_radius - self.center_disk_inner_radius);

        let u = 0.5 + 0.5 * r * frame.cos_phi;
        let v = 0.5 + 0.5 * r * frame.sin_phi;
        UVCoordinates { u, v }
    }

    fn raymarch(
        &self,
        geometry: &dyn Geometry,
        frequency: &RayFrequencyData,
        remaining_steps: &[Step],
    ) -> Result<CIETristimulus, RaytracerError> {
        let first = remaining_steps
            .first()
            .expect("raymarch requires the fired step window");
        trace!("Start raymarching at {:?} ", first);
        let mut segment_state = SegmentState {
            distance_accumulated: 0.0,
            transparency: 1.0,
            alpha_weighted_sum: 0.0,
            alpha_weight_total: 0.0,
            color_accumulated: CIETristimulus::new(0.0, 0.0, 0.0, 0.0),
        };

        // If the first point is outside of the disc it is an entering ray, otherwise it is an exiting ray.
        let is_entering = !self.is_inside_bounding_cylinder(&first.spatial_cartesian());
        if !is_entering {
            // TODO: a camera INSIDE the gas also starts inside the cylinder
            // and is suppressed here; supporting that case needs the window
            // index (entry-guard exception for the ray's first window).
            trace!("Ray is exiting the disc, skipping raymarching.");
            return Ok(CIETristimulus::new(0.0, 0.0, 0.0, 0.0));
        }

        // One call owns exactly one gas episode: it ends at the first
        // inside-to-outside transition of the bounding cylinder, not only at
        // a fully non-overlapping window. Otherwise a hole narrower than the
        // local window spacing (every window starts inside or crosses a
        // wall) would let this call march the far-side episode that the
        // scene's re-entry firing marches again (double counting).
        let mut entered_gas = false;
        for steps in remaining_steps.windows(2) {
            let segment = Segment::from_steps(&steps[0], &steps[1]);
            let from = &segment.start_cartesian;
            let to = &segment.end_cartesian;

            let crosses_boundary = self.intersects_internal(&segment).is_some();
            let overlaps = self.is_inside_bounding_cylinder(from) || crosses_boundary;

            // If the segment does not intersect the disc, we can skip it entirely.
            if !overlaps {
                break;
            }

            let result =
                self.raymarch_segment(from, to, geometry, frequency, &mut segment_state)?;

            if result == RaymarchResult::AlreadyOpaque {
                break;
            }

            let end_inside = self.is_inside_bounding_cylinder(to);
            if !end_inside && (entered_gas || crosses_boundary) {
                // The ray has left the volume: episode over. A later
                // re-entry fires its own call (entry guard lets it through,
                // since it starts outside).
                break;
            }
            entered_gas |= end_inside;
        }

        // Final alpha combines physical opacity with texture alpha, applied once at the end.
        let physical_opacity = 1.0 - segment_state.transparency;
        let texture_alpha = if segment_state.alpha_weight_total > 0.0 {
            segment_state.alpha_weighted_sum / segment_state.alpha_weight_total
        } else {
            1.0
        };
        segment_state.color_accumulated.alpha = physical_opacity * texture_alpha;

        trace!("  resulting color: {:?}", segment_state.color_accumulated);
        Ok(segment_state.color_accumulated)
    }

    fn march_constant_step(
        &self,
        p: &Vector3<f64>,
        geometry: &dyn Geometry,
        frequency: &RayFrequencyData,
        segment_state: &mut SegmentState,
    ) -> Result<(), RaytracerError> {
        let sigma_a = self.absorption;
        let sigma_s = self.scattering;
        let step_size = self.step_size;
        // Resolve the sample into the disc's cylindrical frame once and reuse
        // it for the density, the temperature radius and the UV lookup.
        let frame = self.cylindrical_frame(p);
        let density = self.compute_density_in_frame(&frame);

        if density > 0.0 {
            let sigma_t = sigma_a + sigma_s;
            let tau_cell = step_size * density * sigma_t;
            let cell_transmittance = (-tau_cell).exp();

            // Per-sample redshift from the ray's conserved (p_t, p_phi)
            // and the local circular-orbit Killing coefficients:
            // u.p = u^t p_t + u^phi p_phi, exact at every sample with no
            // parallel transport (see docs/plan-01-per-sample-redshift.md).
            // Where no timelike circular orbit exists the gas is
            // unphysical anyway: it still attenuates (below) but emits
            // nothing.
            let sample_position = Point::new_cartesian(0.0, p[0], p[1], p[2]);
            if let Ok(coefficients) = geometry.circular_orbit_killing_coefficients(&sample_position)
            {
                let emitter_energy =
                    coefficients.u_t * frequency.p_t + coefficients.u_phi * frequency.p_phi;
                let redshift = frequency.observer_energy / emitter_energy;

                let temperature = self.temperature_computer.compute_temperature(frame.r)?;
                let uv = self.get_uv_in_frame(&frame);
                let light_color = self.texture_mapper.color_at_uv(
                    &uv,
                    &crate::rendering::texture::TemperatureData {
                        temperature,
                        redshift,
                    },
                )?;

                // Stefan-Boltzmann law: emission intensity scales with T^4.
                // Use a reference temperature for normalization to boost brightness.
                let intensity_factor =
                    (temperature / self.brightness_reference_temperature).powi(4);

                // Kirchhoff: thermal emissivity couples to absorption,
                // j = sigma_a * rho * B(T). With scattering present but
                // no external illumination the source function is the
                // albedo-weighted Planck function S = (sigma_a/sigma_t) B;
                // integrating a constant source across the cell gives the
                // exact per-cell weight transparency * (S/B) * (1 - e^-tau)
                // (~ sigma_a * rho * ds for thin cells), evaluated BEFORE
                // this cell's own transmittance is applied so the cell
                // does not absorb its own emission twice.
                let emission_weight = if sigma_t > 0.0 {
                    segment_state.transparency * (sigma_a / sigma_t) * (1.0 - cell_transmittance)
                } else {
                    0.0
                };
                let step_emission = light_color.mul_color_part(emission_weight * intensity_factor);

                // Keep texture alpha influence separate from per-step emission accumulation.
                let alpha_sample_weight = density * step_size;
                segment_state.alpha_weighted_sum +=
                    light_color.alpha.clamp(0.0, 1.0) * alpha_sample_weight;
                segment_state.alpha_weight_total += alpha_sample_weight;

                segment_state.color_accumulated.x += step_emission.x;
                segment_state.color_accumulated.y += step_emission.y;
                segment_state.color_accumulated.z += step_emission.z;
            } else {
                // No timelike circular orbit here: unphysical gas; it
                // attenuates (below) but emits nothing.
                trace!("  no timelike circular orbit at {:?}; emission skipped", p);
            }

            segment_state.transparency *= cell_transmittance;
        }
        Ok(())
    }

    // Constant-step marching with per-cell Beer-Lambert transmittance; the
    // basic pattern follows the classic tutorial treatment
    // (https://www.scratchapixel.com/lessons/3d-basic-rendering/volume-rendering-for-developers/ray-marching-get-it-right.html),
    // but the emission is the thermal GRRT source term (see the Kirchhoff
    // comment in march_constant_step and the references in
    // docs/plan-02-segment-marching.md), and the march walks the ray's
    // geodesic step windows: step sizes do not align with the segment
    // boundaries, so the sampling phase (distance_accumulated) carries
    // across segments to keep one uniform comb along the whole episode.
    fn raymarch_segment(
        &self,
        from: &Vector3<f64>,
        to: &Vector3<f64>,
        geometry: &dyn Geometry,
        frequency: &RayFrequencyData,
        segment_state: &mut SegmentState,
    ) -> Result<RaymarchResult, RaytracerError> {
        let segment = to - from;
        let segment_length = segment.norm();
        let direction = segment / segment_length;

        // Start the distance at the accumulated distance from the previous segment.
        let mut dist = segment_state.distance_accumulated;

        while dist < segment_length {
            let p = from + direction * dist;

            self.march_constant_step(&p, geometry, frequency, segment_state)?;

            if segment_state.transparency <= TRANSPARENCY_EARLY_EXIT {
                return Ok(RaymarchResult::AlreadyOpaque);
            }

            dist += self.step_size;
        }

        // Store the remaining distance that was not processed in this segment for the next segment.
        segment_state.distance_accumulated = dist - segment_length;

        Ok(RaymarchResult::Continue)
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

        // Fixed-size accumulator: at most two roots, and this runs on every
        // step window of every ray, so a heap `Vec` here is pure allocator
        // traffic.
        let mut hits = [0.0f64; 2];
        let mut hit_count = 0usize;
        for dist in [t1_dist, t2_dist] {
            let t = dist / segment_length;
            if (0.0..=1.0).contains(&t) {
                let p = from + t * (to - from);
                if p.dot(axis).abs() <= half_height {
                    hits[hit_count] = t;
                    hit_count += 1;
                }
            }
        }

        match hit_count {
            0 => NoIntersection,
            1 => OneIntersection(hits[0]),
            _ => TwoIntersections(hits[0].min(hits[1]), hits[0].max(hits[1])),
        }
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
        // At most six candidate parameters (two tubes x two roots, two caps).
        // Stack-allocated for the same reason as in the clipped-cylinder
        // solver above: this is per-step, per-ray code.
        let mut hits = [0.0f64; 6];
        let mut hit_count = 0usize;
        fn push(t: f64, hits: &mut [f64; 6], hit_count: &mut usize) {
            hits[*hit_count] = t;
            *hit_count += 1;
        }
        let dir = to - from;

        // Use 3.0 * thickness to capture the Gaussian tail
        let capture_height = self.thickness * CAPTURE_HEIGHT_FACTOR;

        // Check Outer Tube
        match self.intersects_clipped_cylinder(from, to, outer_radius, axis, capture_height) {
            OneIntersection(t) => push(t, &mut hits, &mut hit_count),
            TwoIntersections(t1, t2) => {
                push(t1, &mut hits, &mut hit_count);
                push(t2, &mut hits, &mut hit_count);
            }
            _ => {}
        }
        // Check Inner Tube
        match self.intersects_clipped_cylinder(from, to, inner_radius, axis, capture_height) {
            OneIntersection(t) => push(t, &mut hits, &mut hit_count),
            TwoIntersections(t1, t2) => {
                push(t1, &mut hits, &mut hit_count);
                push(t2, &mut hits, &mut hit_count);
            }
            _ => {}
        }
        // Check Caps
        for pos in [capture_height, -capture_height] {
            if let Some(t) = self.intersects_cap(from, to, outer_radius, axis, pos) {
                let p = from + t * dir;
                let r_sq = p.cross(axis).norm_squared();
                if r_sq >= inner_radius * inner_radius {
                    push(t, &mut hits, &mut hit_count);
                }
            }
        }

        let hits = &mut hits[..hit_count];
        hits.sort_by(|a, b| a.total_cmp(b));

        match hit_count {
            0 => NoIntersection,
            1 => OneIntersection(hits[0]),
            _ => TwoIntersections(hits[0], hits[1]),
        }
    }

    fn is_inside_bounding_cylinder(&self, p: &Vector3<f64>) -> bool {
        if p.dot(&self.axis).abs() > self.thickness * CAPTURE_HEIGHT_FACTOR {
            return false;
        }
        // Compared squared to skip the square root; `axis` is a unit vector,
        // so |p x axis|^2 = |p|^2 - (p . axis)^2.
        let axial = p.dot(&self.axis);
        let radial_squared = (p.norm_squared() - axial * axial).max(0.0);
        radial_squared > self.center_disk_inner_radius * self.center_disk_inner_radius
            && radial_squared < self.center_disk_outer_radius * self.center_disk_outer_radius
    }
}

#[derive(Debug)]
pub enum CylinderIntersection {
    NoIntersection,
    Parallel,
    OneIntersection(f64),
    TwoIntersections(f64, f64),
}

impl VolumetricDisc {
    // The volume boundary is deliberately a Cartesian cylinder (the disc is a
    // stylized visual model); it does not need the metric radial coordinate.
    fn intersects_internal(&self, segment: &Segment) -> Option<Intersection> {
        let y_start_spatial = segment.start_cartesian;
        let y_end_spatial = segment.end_cartesian;

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
        })
    }
}

impl Hittable for VolumetricDisc {
    fn intersects(&self, segment: &Segment, _geometry: &dyn Geometry) -> Option<Intersection> {
        self.intersects_internal(segment)
    }

    fn bounding_radius(&self, _geometry: &dyn Geometry) -> f64 {
        // The volume boundary is a Cartesian cylinder (see
        // `intersects_internal`), so its corner is the farthest reachable
        // point.
        let capture_height = self.thickness * CAPTURE_HEIGHT_FACTOR;
        (self.center_disk_outer_radius * self.center_disk_outer_radius
            + capture_height * capture_height)
            .sqrt()
    }

    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData,
        geometry: &dyn Geometry,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.raymarch(
            geometry,
            &color_computation_data.frequency,
            color_computation_data.remaining_steps,
        )
    }

    // The raymarch computes redshift and temperature PER SAMPLE from the
    // conserved ray quantities; the per-intersection values built from these
    // two queries are never read for volumetric hits. Returning fixed values
    // also removes a spurious per-pixel failure path: the real queries error
    // near the hole (NoCircularOrbitPossible below the photon orbit,
    // BelowRISCO from the temperature LUT) for entry crossings at the inner
    // wall, aborting whole pixels for values nobody consumes.
    fn energy_of_emitter(
        &self,
        _geometry: &dyn Geometry,
        _step: &Step,
    ) -> Result<f64, RaytracerError> {
        Ok(1.0)
    }

    fn temperature_of_emitter(
        &self,
        _point: &Point,
        _geometry: &dyn Geometry,
    ) -> Result<f64, RaytracerError> {
        Ok(1.0)
    }
}

impl SceneObject for VolumetricDisc {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::rendering::color::Color;
    use crate::rendering::texture::{CheckerMapper, TemperatureData, TextureMap};
    use approx::assert_abs_diff_eq;
    use std::sync::Arc;

    fn unit_frequency() -> RayFrequencyData {
        RayFrequencyData {
            observer_energy: 1.0,
            p_t: 1.0,
            p_phi: 0.0,
        }
    }

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
        ) -> Result<CIETristimulus, RaytracerError> {
            Ok(self.color)
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

        assert!(
            disc.intersects_internal(&Segment::from_points(&y_start, &y_end))
                .is_some()
        );
    }

    #[test]
    fn test_volumetric_disc_intersection_miss_above_caps() {
        let disc = create_disc();
        let y_start = Point::new_cartesian(0.0, 1.5, 0.0, 2.0);
        let y_end = Point::new_cartesian(0.0, 2.5, 0.0, 2.0);

        assert!(
            disc.intersects_internal(&Segment::from_points(&y_start, &y_end))
                .is_none()
        );
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

    /// Two-point straight step slice mimicking a fired entry window: starts
    /// outside the bounding cylinder (inner hole) and crosses the gas.
    fn straight_steps(from: Vector3<f64>, to: Vector3<f64>) -> Vec<Step> {
        let mk = |v: &Vector3<f64>| {
            Step::new(
                Point::new_cartesian(0.0, v[0], v[1], v[2]),
                FourVector::new_cartesian(1.0, 1.0, 0.0, 0.0),
                0.0,
                0,
            )
        };
        vec![mk(&from), mk(&to)]
    }

    #[test]
    fn test_volumetric_disc_raymarch_produces_opacity() {
        let disc = create_disc();
        let steps = straight_steps(Vector3::new(0.2, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        let color = disc
            .raymarch(&EuclideanSpace::new(), &unit_frequency(), &steps)
            .expect("raymarch should succeed");

        assert!(color.alpha > 0.0);
        assert!(color.alpha <= 1.0);
    }

    /// Kirchhoff's law: a purely scattering medium (sigma_a = 0) redirects
    /// light but cannot create thermal photons; it must attenuate (alpha > 0)
    /// while emitting nothing.
    #[test]
    fn test_pure_scattering_gas_attenuates_but_does_not_emit() {
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
            0.0, // absorption: none
            0.4, // scattering only
            Vector3::new(1.0, 1.0, 1.0),
            1.0,
        );
        let steps = straight_steps(Vector3::new(0.2, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        let color = disc
            .raymarch(&EuclideanSpace::new(), &unit_frequency(), &steps)
            .expect("raymarch should succeed");

        assert_eq!(color.x, 0.0);
        assert_eq!(color.y, 0.0);
        assert_eq!(color.z, 0.0);
        assert!(color.alpha > 0.0);
    }

    /// Fully transparent coefficients (sigma_a = sigma_s = 0) must yield no
    /// interaction at all, not NaN from 0/0 in the source-function weight.
    #[test]
    fn test_zero_coefficients_yield_no_interaction() {
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
            0.0,
            0.0,
            Vector3::new(1.0, 1.0, 1.0),
            1.0,
        );
        let steps = straight_steps(Vector3::new(0.2, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        let color = disc
            .raymarch(&EuclideanSpace::new(), &unit_frequency(), &steps)
            .expect("raymarch should succeed");

        assert!(color.x.is_finite() && color.x == 0.0);
        assert_eq!(color.alpha, 0.0);
    }

    /// Segmentation invariance: marching the same straight path as ONE
    /// window must equal marching it split into MANY windows with
    /// boundaries deliberately misaligned to the step size. This pins the
    /// phase carry-over (distance_accumulated): if a segment boundary ever
    /// reset or double-counted the sampling comb, the results diverge.
    #[test]
    fn test_raymarch_segmentation_invariance() {
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

        let one = straight_steps(Vector3::new(0.2, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        // Same path cut at odd places; the first window still contains the
        // entry crossing (as a fired window would), later cuts fall inside
        // the gas and one straddles the outer wall.
        let xs = [0.2, 1.35, 2.1, 2.8005, 3.4, 5.0];
        let many: Vec<Step> = xs
            .iter()
            .map(|x| {
                Step::new(
                    Point::new_cartesian(0.0, *x, 0.0, 0.0),
                    FourVector::new_cartesian(1.0, 1.0, 0.0, 0.0),
                    0.0,
                    0,
                )
            })
            .collect();

        let geometry = EuclideanSpace::new();
        let a = disc
            .raymarch(&geometry, &unit_frequency(), &one)
            .expect("single-window raymarch should succeed");
        let b = disc
            .raymarch(&geometry, &unit_frequency(), &many)
            .expect("multi-window raymarch should succeed");

        assert_abs_diff_eq!(a.x, b.x, epsilon = 1e-9);
        assert_abs_diff_eq!(a.y, b.y, epsilon = 1e-9);
        assert_abs_diff_eq!(a.z, b.z, epsilon = 1e-9);
        assert_abs_diff_eq!(a.alpha, b.alpha, epsilon = 1e-9);
    }

    /// Entry guard: a call whose slice starts INSIDE the bounding cylinder
    /// is an exit/interior firing; a previous entry call already owns that
    /// episode, so this call must contribute nothing (double-count guard).
    #[test]
    fn test_exit_firing_is_suppressed() {
        let disc = create_disc();
        // (2,0,0) lies inside the gas annulus (1 < r < 3, h = 0).
        let steps = straight_steps(Vector3::new(2.0, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        let color = disc
            .raymarch(&EuclideanSpace::new(), &unit_frequency(), &steps)
            .expect("raymarch should succeed");

        assert_eq!(color.x, 0.0);
        assert_eq!(color.y, 0.0);
        assert_eq!(color.z, 0.0);
        assert_eq!(color.alpha, 0.0);
    }

    /// Episode protocol: a chord through the inner hole crosses the gas
    /// twice. One call owns exactly the FIRST episode (the driver breaks in
    /// the hole), and the re-entry is handled by a fresh call on the tail,
    /// mirroring the scene's re-firing on the re-entry crossing.
    #[test]
    fn test_one_episode_per_call_with_reentry() {
        let disc = create_disc();
        let mk = |x: f64| {
            Step::new(
                Point::new_cartesian(0.0, x, 0.5, 0.0),
                FourVector::new_cartesian(1.0, 1.0, 0.0, 0.0),
                0.0,
                0,
            )
        };
        // Chord at y = 0.5: gas episodes at |x| in ~(0.866, 2.958); the
        // window (-0.5 -> 0.5) lies wholly in the inner hole.
        let full: Vec<Step> = [-5.0, -2.0, -0.5, 0.5, 2.0, 5.0]
            .map(mk)
            .into_iter()
            .collect();
        let episode1: Vec<Step> = [-5.0, -2.0, -0.5, 0.5].map(mk).into_iter().collect();
        let episode2: Vec<Step> = [0.5, 2.0, 5.0].map(mk).into_iter().collect();

        let geometry = EuclideanSpace::new();
        let a = disc
            .raymarch(&geometry, &unit_frequency(), &full)
            .expect("full-slice raymarch should succeed");
        let b = disc
            .raymarch(&geometry, &unit_frequency(), &episode1)
            .expect("episode-1 raymarch should succeed");
        let c = disc
            .raymarch(&geometry, &unit_frequency(), &episode2)
            .expect("episode-2 raymarch should succeed");

        // The full slice accumulates ONLY episode 1 (break in the hole).
        assert_eq!(a.x, b.x);
        assert_eq!(a.y, b.y);
        assert_eq!(a.alpha, b.alpha);
        // The re-entry call (starting outside, in the hole) marches episode 2.
        assert!(c.alpha > 0.0);
    }

    /// PR #121 review (Codex): if the inner hole is narrower than the local
    /// window spacing, NO window lies wholly outside the volume - each one
    /// starts in gas or crosses a wall. The call must still end at the first
    /// inside-to-outside transition; otherwise it marches the far-side
    /// episode that the re-entry firing marches again (double counting).
    #[test]
    fn test_narrow_hole_does_not_double_count_far_episode() {
        let disc = create_disc();
        let mk = |x: f64| {
            Step::new(
                Point::new_cartesian(0.0, x, 0.5, 0.0),
                FourVector::new_cartesian(1.0, 1.0, 0.0, 0.0),
                0.0,
                0,
            )
        };
        // Chord at y = 0.5: gas at |x| in ~(0.866, 2.958). The cut at -0.7
        // (r ~ 0.86) is the only point in the hole; the window (-0.7 -> 0.9)
        // crosses back INTO gas, so no window is fully outside.
        let full: Vec<Step> = [-5.0, -2.0, -0.7, 0.9, 2.0, 5.0]
            .map(mk)
            .into_iter()
            .collect();
        let episode1: Vec<Step> = [-5.0, -2.0, -0.7].map(mk).into_iter().collect();
        let reentry: Vec<Step> = [-0.7, 0.9, 2.0, 5.0].map(mk).into_iter().collect();

        let geometry = EuclideanSpace::new();
        let a = disc
            .raymarch(&geometry, &unit_frequency(), &full)
            .expect("full-slice raymarch should succeed");
        let b = disc
            .raymarch(&geometry, &unit_frequency(), &episode1)
            .expect("episode-1 raymarch should succeed");
        let c = disc
            .raymarch(&geometry, &unit_frequency(), &reentry)
            .expect("re-entry raymarch should succeed");

        // The full-slice call ends at the gas -> hole transition; the far
        // episode belongs exclusively to the re-entry call.
        assert_eq!(a.x, b.x);
        assert_eq!(a.alpha, b.alpha);
        assert!(c.alpha > 0.0);
    }

    /// Optically thick gas saturates: opacity approaches 1 (and the marcher
    /// is allowed to stop early once it does).
    #[test]
    fn test_dense_gas_saturates_opacity() {
        let disc = VolumetricDisc::new(
            1.0,
            3.0,
            Arc::new(FixedTextureMap {
                color: CIETristimulus::new(2.0, 1.0, 0.5, 1.0),
            }),
            Box::new(DummyTemperatureComputer),
            Vector3::new(0.0, 0.0, 1.0),
            4,
            42,
            500,
            0.01,
            0.5,
            1e6, // very dense
            1000.0,
            5.0,
            5.0,
            Vector3::new(1.0, 1.0, 1.0),
            1.0,
        );
        let steps = straight_steps(Vector3::new(0.2, 0.0, 0.0), Vector3::new(5.0, 0.0, 0.0));

        let color = disc
            .raymarch(&EuclideanSpace::new(), &unit_frequency(), &steps)
            .expect("raymarch should succeed");

        assert!(color.alpha > 0.99, "alpha = {}", color.alpha);
    }
}
