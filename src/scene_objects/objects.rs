use crate::geometry::coordinates::conversion::CartesianConvertible;
use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::Step;
use crate::rendering::light::radiance::Radiance;
use crate::rendering::light::redshift::{RayFrequencyData, RedshiftComputer};
use crate::rendering::texture::TemperatureData;
use crate::scene_objects::hittable::{ColorComputationData, Hittable};

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
fn step_at_intersection<C: CartesianConvertible>(
    y_start: &Step<C>,
    y_end: &Step<C>,
    intersection_point: Point<C>,
    t: f64,
) -> Step<C> {
    // Ensure all the chart enums align.
    debug_assert_eq!(y_start.x.chart, y_end.x.chart);
    debug_assert_eq!(y_start.p.at.chart, y_end.p.at.chart);
    debug_assert_eq!(y_start.x.chart, intersection_point.chart);
    let s = 1.0 - t;
    Step {
        t: s * y_start.t + t * y_end.t,
        step: y_start.step,
        x: intersection_point,
        p: FourVector::new(
            s * y_start.p[0] + t * y_end.p[0],
            s * y_start.p[1] + t * y_end.p[1],
            s * y_start.p[2] + t * y_end.p[2],
            s * y_start.p[3] + t * y_end.p[3],
            intersection_point,
        ),
    }
}

pub struct Objects<'a, G: Geometry> {
    geometry: &'a G,
    objects: Vec<Box<dyn Hittable<G::Chart>>>,
}

impl<'a, G: Geometry> Objects<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self {
            geometry,
            objects: Vec::new(),
        }
    }

    pub fn add_object(&mut self, hittable: Box<dyn Hittable<G::Chart>>) {
        self.objects.push(hittable);
    }

    pub fn intersects(
        &self,
        y_start: &Step<G::Chart>,
        y_end: &Step<G::Chart>,
        frequency: &RayFrequencyData,
        remaining_steps: &[Step<G::Chart>],
    ) -> Result<Option<Radiance>, RaytracerError> {
        // Right after a chart switch the two ends of a segment are in different charts. Express
        // the start in the end's chart, so that objects interpolate within one chart.
        let aligned_start;
        let y_start = if y_start.x.chart != y_end.x.chart {
            aligned_start = y_start.in_chart(y_end.x.chart);
            &aligned_start
        } else {
            y_start
        };

        let redshift_computer = RedshiftComputer::new(self.geometry);
        let mut resulting_color = None;
        let mut shortest_distance = f64::MAX;

        let y_start_cartesian = y_start.x.get_spatial_vector_cartesian();

        // The step size can be rather large, so it makes sense to sort the objects by their
        // distance to the y_start point.
        for hittable in &self.objects {
            if let Some(intersection_data) = hittable.intersects(&y_start, &y_end, self.geometry) {
                let intersection_point = intersection_data
                    .intersection_point
                    .get_spatial_vector_cartesian();
                let distance = (intersection_point - y_start_cartesian).norm();
                if distance < shortest_distance {
                    shortest_distance = distance;
                    let intersection_point = G::Chart::point_from_cartesian(
                        y_start.x.chart,
                        &intersection_data.intersection_point,
                    );
                    let intersection_step = step_at_intersection(
                        y_start,
                        y_end,
                        intersection_point,
                        intersection_data.t,
                    );
                    let emitter_energy =
                        hittable.energy_of_emitter(self.geometry, &intersection_step)?;
                    let redshift = redshift_computer
                        .compute_redshift_from_energies(emitter_energy, frequency.observer_energy);
                    let temperature =
                        hittable.temperature_of_emitter(&intersection_point, self.geometry)?;

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
    use crate::geometry::coordinates::cartesian::CartesianChart;
    use crate::geometry::spacetimes::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::point::Point;
    use crate::rendering::light::color::Color;
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

    fn create_sphere_at(
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        color: u8,
    ) -> Box<dyn Hittable<CartesianChart>> {
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
    fn test_add_and_intersect_two_sphere() {
        let geometry = EuclideanSpace::new();
        let mut objects = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 1.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects.add_object(farther_sphere);
        objects.add_object(closer_sphere);

        let step_start = Step {
            t: 0.0,
            step: 0,
            x: Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            p: FourVector::new(
                1.0,
                0.0,
                0.0,
                0.0,
                Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            ),
        };
        let step_end = Step {
            t: 100.0,
            step: 1000,
            x: Point::new_cartesian(0.0, 0.0, 0.0, 3.0),
            p: FourVector::new(1.0, 0.0, 0.0, 0.0, Point::new_cartesian(0.0, 0.0, 0.0, 3.0)),
        };

        let result = objects
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_add_and_intersect_spheres_inside_each_other() {
        let step_start = Step {
            t: 0.0,
            step: 0,
            x: Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            p: FourVector::new(
                1.0,
                0.0,
                0.0,
                0.0,
                Point::new_cartesian(0.0, 0.0, 0.0, -3.0),
            ),
        };
        let step_end = Step {
            t: 100.0,
            step: 1000,
            x: Point::new_cartesian(0.0, 0.0, 0.0, 0.0),
            p: FourVector::new(1.0, 0.0, 0.0, 0.0, Point::new_cartesian(0.0, 0.0, 0.0, 0.0)),
        };

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

    #[test]
    fn test_intersect_disc_with_schwarzschild_native_spherical_steps() {
        // Regression test: Schwarzschild's native coordinate system is
        // Spherical, but Disc::intersects always tags its returned
        // intersection point Cartesian. An inclined (non-equatorial) ray
        // that crosses the disc plane must not panic or silently misread
        // the mismatched coordinate tags when building the emitter step.
        use crate::geometry::spacetimes::schwarzschild::Schwarzschild;
        use crate::rendering::light::temperature::ConstantTemperatureComputer;
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
            1.0,
        )));

        // r = 6 (inside the disc annulus), phi = 0, theta straddling the
        // equatorial plane pi/2 where the disc lives.
        let step_start = Step {
            t: 0.0,
            step: 0,
            x: Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 - 0.3, 0.0),
            p: FourVector::new(
                1.0,
                0.0,
                0.1,
                0.0,
                Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 - 0.3, 0.0),
            ),
        };
        let step_end = Step {
            t: 1.0,
            step: 1,
            x: Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 + 0.3, 0.0),
            p: FourVector::new(
                1.0,
                0.0,
                0.1,
                0.0,
                Point::new_spherical(0.0, 6.0, std::f64::consts::FRAC_PI_2 + 0.3, 0.0),
            ),
        };

        let result = objects
            .intersects(&step_start, &step_end, &unit_frequency(), &[])
            .unwrap();
        assert!(result.is_some());

        // The same segment right after a chart switch: the end in the rotated chart.
        use crate::geometry::coordinates::spherical::RotatedSpherical;
        use crate::geometry::atlases::spherical::SphericalChart;
        let rotated = SphericalChart::Rotated(RotatedSpherical { rotation: 90.0 });
        let step_end_rotated = step_end.in_chart(rotated);
        assert_ne!(step_start.x.chart, step_end_rotated.x.chart);
        let switched = objects
            .intersects(&step_start, &step_end_rotated, &unit_frequency(), &[])
            .unwrap();
        let (expected, switched) = (result.unwrap(), switched.unwrap());
        assert_abs_diff_eq!(switched.x, expected.x, epsilon = 1e-9);
        assert_abs_diff_eq!(switched.y, expected.y, epsilon = 1e-9);
        assert_abs_diff_eq!(switched.z, expected.z, epsilon = 1e-9);
    }
}
