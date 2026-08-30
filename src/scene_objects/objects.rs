use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::redshift::{RayFrequencyData, RedshiftComputer};
use crate::rendering::texture::TemperatureData;
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Segment};
use nalgebra::Vector3;

/// Build the integration `Step` at the exact object intersection: the position
/// is the already-solved intersection point (exact, not interpolated), and the
/// momentum is linearly interpolated between `y_start` (t = 0) and `y_end`
/// (t = 1). The momentum has no closed-form solve at the intersection, so
/// linear interpolation in the integrator's native coordinate system is the
/// best available estimate; using the exact intersection point for position
/// avoids the mismatch that comes from interpolating curvilinear coordinates
/// (r, theta, phi) along what is really a straight line in Cartesian space.
///
/// `Hittable::intersects` implementations are free to tag their returned
/// point in whatever convention is easiest for that object (e.g. `Sphere`
/// always returns `Spherical`, `Disc` always returns `Cartesian`), which does
/// not necessarily match the geometry's own native coordinate system that
/// `energy_of_emitter`'s inner products and velocity fields assume. Convert
/// the intersection point into `y_start`/`y_end`'s coordinate system (the
/// geometry's native one) before using it as the step's position.
fn step_at_intersection(y_start: &Step, y_end: &Step, intersection_point: Point, t: f64) -> Step {
    debug_assert_eq!(y_start.x.coordinate_system, y_end.x.coordinate_system);
    debug_assert_eq!(y_start.p.coordinate_system, y_end.p.coordinate_system);
    let cs = y_start.p.coordinate_system;
    let s = 1.0 - t;
    Step::new(
        intersection_point.to_coordinate_system(y_start.x.coordinate_system),
        FourVector::new(
            s * y_start.p[0] + t * y_end.p[0],
            s * y_start.p[1] + t * y_end.p[1],
            s * y_start.p[2] + t * y_end.p[2],
            s * y_start.p[3] + t * y_end.p[3],
            cs,
        ),
        s * y_start.t + t * y_end.t,
        y_start.step,
    )
}

pub trait SceneObject: Hittable {}

pub struct Objects<'a, G: Geometry> {
    geometry: &'a G,
    objects: Vec<Box<dyn SceneObject>>,
    redshift_computer: RedshiftComputer<'a, G>,
    /// Squared radius of an origin-centred sphere containing every object,
    /// maintained as objects are added. See `chord_reaches_objects`.
    bounding_radius_squared: f64,
}

/// Squared distance from the origin to the nearest point of the straight
/// segment `start` -> `end`.
fn closest_approach_squared(start: &Vector3<f64>, end: &Vector3<f64>) -> f64 {
    let direction = end - start;
    let length_squared = direction.norm_squared();
    if length_squared <= 0.0 {
        return start.norm_squared();
    }
    let t = (-start.dot(&direction) / length_squared).clamp(0.0, 1.0);
    (start + t * direction).norm_squared()
}

impl<'a, G: Geometry> Objects<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self {
            geometry,
            objects: Vec::new(),
            redshift_computer: RedshiftComputer::new(geometry),
            bounding_radius_squared: 0.0,
        }
    }

    pub fn add_object(&mut self, hittable: Box<dyn SceneObject>) {
        let radius = hittable.bounding_radius(self.geometry);
        self.bounding_radius_squared = self.bounding_radius_squared.max(radius * radius);
        self.objects.push(hittable);
    }

    /// Broad phase: can this step window's chord reach any object at all?
    ///
    /// Every object tests the straight chord between the two steps, so if the
    /// chord never comes within the scene's bounding radius of the origin, no
    /// object can report a hit on it. Rays that escape to the celestial sphere
    /// spend nearly all of their thousands of steps out there, and this turns
    /// each of those steps from one narrow-phase test per object into a
    /// handful of arithmetic.
    fn chord_reaches_objects(&self, segment: &Segment) -> bool {
        closest_approach_squared(&segment.start_cartesian, &segment.end_cartesian)
            <= self.bounding_radius_squared
    }

    pub fn intersects(
        &self,
        y_start: &Step,
        y_end: &Step,
        frequency: &RayFrequencyData,
        remaining_steps: &[Step],
    ) -> Result<Option<CIETristimulus>, RaytracerError> {
        let mut resulting_color = None;
        let mut shortest_distance_squared = f64::MAX;

        let segment = Segment::from_steps(y_start, y_end);
        if !self.chord_reaches_objects(&segment) {
            return Ok(None);
        }
        let y_start_cartesian = segment.start_cartesian;

        // The step size can be rather large, so it makes sense to sort the objects by their
        // distance to the y_start point.
        for hittable in &self.objects {
            if let Some(intersection_data) = hittable.intersects(&segment, self.geometry) {
                let intersection_point = intersection_data
                    .intersection_point
                    .get_spatial_vector_cartesian();
                // Compared squared: the ordering is the same and the square
                // root is pure cost in a per-step loop.
                let distance_squared = (intersection_point - y_start_cartesian).norm_squared();
                if distance_squared < shortest_distance_squared {
                    shortest_distance_squared = distance_squared;
                    let intersection_step = step_at_intersection(
                        y_start,
                        y_end,
                        intersection_data.intersection_point,
                        intersection_data.t,
                    );
                    let emitter_energy =
                        hittable.energy_of_emitter(self.geometry, &intersection_step)?;
                    let redshift = self
                        .redshift_computer
                        .compute_redshift_from_energies(emitter_energy, frequency.observer_energy);
                    let temperature = hittable.temperature_of_emitter(
                        &intersection_data.intersection_point,
                        self.geometry,
                    )?;

                    let color_computation_data = ColorComputationData {
                        uv: intersection_data.uv,
                        temperature_data: TemperatureData {
                            redshift,
                            temperature,
                        },
                        frequency: *frequency,
                        remaining_steps,
                    };
                    resulting_color =
                        Some(hittable.color_at_uv(&color_computation_data, self.geometry)?);
                }
            }
        }
        Ok(resulting_color)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::color::Color;
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects::sphere::Sphere;
    use approx::assert_abs_diff_eq;
    use std::sync::Arc;

    /// Flat-space stand-in for the per-ray conserved quantities: an emitter
    /// at rest sees g = observer_energy / p_t = 1.
    fn unit_frequency() -> RayFrequencyData {
        RayFrequencyData {
            observer_energy: 1.0,
            p_t: 1.0,
            p_phi: 0.0,
        }
    }

    fn create_sphere_at(x: f64, y: f64, z: f64, radius: f64, color: u8) -> Box<dyn SceneObject> {
        Box::new(Sphere::new(
            radius,
            Arc::new(CheckerMapper::new(
                3.0,
                5.0,
                5.0,
                Color::new(color, color, color, 255),
                Color::new(color, color, color, 255),
            )),
            Point::new_cartesian(0.0, x, y, z),
            0.0,
        ))
    }
    #[test]
    #[ignore] // TODO: if step goes through sphere fully, not intersection seen. Need to fix.
    fn test_add_and_intersect_two_sphere() {
        let geometry = EuclideanSpace::new();
        let mut objects = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 1.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects.add_object(farther_sphere);
        objects.add_object(closer_sphere);

        let step_start = Step::new(
            Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            0.0,
            0,
        );
        let step_end = Step::new(
            Point::new_cartesian(0.0, 0.0, 0.0, 3.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            100.0,
            1000,
        );

        let result = objects
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_add_and_intersect_spheres_inside_each_other() {
        let step_start = Step::new(
            Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            0.0,
            0,
        );
        let step_end = Step::new(
            Point::new_cartesian(0.0, 0.0, 0.0, 0.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0),
            100.0,
            1000,
        );

        let geometry = EuclideanSpace::new();
        let mut objects_setup_1 = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects_setup_1.add_object(farther_sphere);
        objects_setup_1.add_object(closer_sphere);

        let result_1 = objects_setup_1
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result_1.is_some());
        assert_abs_diff_eq!(result_1.unwrap().x, 0.121, epsilon = 1e-2);

        let mut objects_setup_2 = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);
        objects_setup_2.add_object(closer_sphere);
        objects_setup_2.add_object(farther_sphere);

        let result_2 = objects_setup_2
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result_2.is_some());
        assert_abs_diff_eq!(result_2.unwrap().x, 0.121, epsilon = 1e-2);
    }

    /// The broad phase must never reject a window that a narrow-phase test
    /// would have accepted. Kerr is the case where a naive bound gets it
    /// wrong: the equatorial relation `x^2 + y^2 = r^2 + a^2` puts the disc's
    /// outer edge at a Cartesian radius LARGER than its (metric) outer
    /// radius, so bounding the scene by the outer radius alone would clip the
    /// rim off every Kerr disc.
    #[test]
    fn test_broad_phase_keeps_kerr_disc_rim_beyond_its_metric_outer_radius() {
        use crate::geometry::kerr::Kerr;
        use crate::rendering::temperature::ConstantTemperatureComputer;
        use crate::scene_objects::disc::Disc;

        let a = 0.499_f64;
        let outer = 8.0_f64;
        let geometry = Kerr::new(1.0, a, 1e-4);
        let disc = Disc::new(
            0.795,
            outer,
            Arc::new(CheckerMapper::new(
                0.0,
                5.0,
                5.0,
                Color::new(100, 100, 100, 255),
                Color::new(200, 200, 200, 255),
            )),
            Box::new(ConstantTemperatureComputer::new(5000.0)),
        );
        let mut objects = Objects::new(&geometry);
        objects.add_object(Box::new(disc));

        // Just inside the metric outer edge, but at a Cartesian radius past
        // `outer` - exactly the band a too-tight bound would discard.
        let x = (outer * outer + a * a).sqrt() - 1e-3;
        assert!(x > outer, "test point must lie outside the metric radius");
        let step_start = Step::new(
            Point::new_cartesian(0.0, x, 0.0, 0.1),
            FourVector::new_cartesian(1.0, 0.0, 0.0, -1.0),
            0.0,
            0,
        );
        let step_end = Step::new(
            Point::new_cartesian(0.0, x, 0.0, -0.1),
            FourVector::new_cartesian(1.0, 0.0, 0.0, -1.0),
            1.0,
            1,
        );

        let segment = Segment::from_steps(&step_start, &step_end);
        assert!(
            objects.chord_reaches_objects(&segment),
            "disc rim at Cartesian radius {x} was rejected by the broad phase"
        );
        // And the narrow phase does find gas there, so the rejection would
        // have lost a real hit.
        assert!(objects.objects[0].intersects(&segment, &geometry).is_some());
    }

    /// ... and it must actually reject: a window that stays far outside the
    /// scene's bound never reaches the narrow phase.
    #[test]
    fn test_broad_phase_rejects_chords_that_stay_outside_the_scene() {
        let geometry = EuclideanSpace::new();
        let mut objects = Objects::new(&geometry);
        objects.add_object(create_sphere_at(0.0, 0.0, 0.0, 1.0, 100));

        // A chord crossing the disc plane, but 50 units away from the sphere.
        let step_start = Step::new(
            Point::new_cartesian(0.0, 50.0, 0.0, 1.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, -1.0),
            0.0,
            0,
        );
        let step_end = Step::new(
            Point::new_cartesian(0.0, 50.0, 0.0, -1.0),
            FourVector::new_cartesian(1.0, 0.0, 0.0, -1.0),
            1.0,
            1,
        );

        assert!(!objects.chord_reaches_objects(&Segment::from_steps(&step_start, &step_end)));
        assert!(
            objects
                .intersects(&step_start, &step_end, &unit_frequency(), &[])
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn test_intersect_disc_with_schwarzschild_native_spherical_steps() {
        // Regression test: Schwarzschild's native coordinate system is
        // Spherical, but Disc::intersects always tags its returned
        // intersection point Cartesian. An inclined (non-equatorial) ray
        // that crosses the disc plane must not panic or silently misread
        // the mismatched coordinate tags when building the emitter step.
        use crate::geometry::schwarzschild::Schwarzschild;
        use crate::rendering::temperature::ConstantTemperatureComputer;
        use crate::scene_objects::disc::Disc;

        let geometry = Schwarzschild::new(2.0, 1e-4);
        let mut objects = Objects::new(&geometry);
        objects.add_object(Box::new(Disc::new(
            4.0,
            10.0,
            Arc::new(CheckerMapper::new(
                3.0,
                5.0,
                5.0,
                Color::new(100, 100, 100, 255),
                Color::new(200, 200, 200, 255),
            )),
            Box::new(ConstantTemperatureComputer::new(5000.0)),
        )));

        // r = 6 (inside the disc annulus), phi = 0, theta straddling the
        // equatorial plane pi/2 where the disc lives.
        let step_start = Step::new(
            Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 - 0.3, 0.0),
            FourVector::new_spherical(1.0, 0.0, 0.1, 0.0),
            0.0,
            0,
        );
        let step_end = Step::new(
            Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 + 0.3, 0.0),
            FourVector::new_spherical(1.0, 0.0, 0.1, 0.0),
            1.0,
            1,
        );

        let result = objects
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result.is_some());
    }
}
