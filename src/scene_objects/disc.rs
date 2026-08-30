use crate::geometry::geometry::Geometry;
use crate::geometry::point::{CoordinateSystem, Point};
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::temperature::TemperatureComputer;
use crate::rendering::texture::{TextureMapHandle, UVCoordinates};
use crate::scene_objects::hittable::{ColorComputationData, Hittable, Intersection};
use crate::scene_objects::objects::SceneObject;
use nalgebra::Vector3;
use roots::{SimpleConvergency, find_root_brent};

pub struct Disc {
    center_disk_inner_radius: f64,
    center_disk_outer_radius: f64,
    texture_mapper: TextureMapHandle,
    temperature_computer: Box<dyn TemperatureComputer>,
}

impl Disc {
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

    fn create_intersection(
        &self,
        geometry: &dyn Geometry,
        center: Vector3<f64>,
        y_start_spatial: Vector3<f64>,
        direction: Vector3<f64>,
        t: f64,
    ) -> Option<Intersection> {
        // In-plane coordinates are interpolated at t; z is pinned to the disc
        // plane (z = 0). For the linear branch the chord already crosses at
        // z = 0, but the cubic branch solves z(s) = 0 on the reconstructed
        // curve while x and y are still taken from the straight chord, so the
        // chord's z at that t is nonzero and must not leak into the radius,
        // UV, emitter velocity or distance.
        let intersection_point = y_start_spatial + t * direction;
        let intersection_point_p =
            Point::new_cartesian(0.0, intersection_point[0], intersection_point[1], 0.0);
        // The disc's inner/outer radius are metric radial coordinates (the
        // convention the ISCO, orbit velocities, and temperature profile use).
        // In Kerr charts the equatorial Cartesian radius obeys
        // x^2 + y^2 = r^2 + a^2, so the Cartesian norm must not be compared
        // against them directly.
        let r = geometry.get_radial_coordinate(&intersection_point_p);

        if r >= self.center_disk_inner_radius && r <= self.center_disk_outer_radius {
            let vector_in_plane = intersection_point - center;

            let phi = vector_in_plane[1].atan2(vector_in_plane[0]); // phi in x-y plane.
            let r_normalized = (r - self.center_disk_inner_radius)
                / (self.center_disk_outer_radius - self.center_disk_inner_radius);

            let u = 0.5 + 0.5 * r_normalized * phi.cos();
            let v = 0.5 + 0.5 * r_normalized * phi.sin();

            return Some(Intersection {
                uv: UVCoordinates { u, v },
                intersection_point: intersection_point_p,
                t,
            });
        } else {
            return None;
        }
    }
}

impl Hittable for Disc {
    fn intersects(
        &self,
        y_start: &Step,
        y_end: &Step,
        geometry: &dyn Geometry,
    ) -> Option<Intersection> {
        // z x y
        let normal = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::new(0.0, 0.0, 0.0);
        let y_start_spatial = y_start.x.get_spatial_vector_cartesian();
        let y_end_spatial = y_end.x.get_spatial_vector_cartesian();
        let direction = y_end_spatial - y_start_spatial;

        // Note: most of the following assumes that normal will be in the z direction.

        // Start and end point along the normal direction.
        let point_start = (center - y_start_spatial).dot(&normal);
        let point_end = direction.dot(&normal);

        // The segment crosses the plane.
        let z0 = y_start_spatial[2];
        let z1 = y_end_spatial[2];
        if z0.signum() != z1.signum() {
            // Opposite z signs mean z0 != z1, so point_end = z1 - z0 is never
            // zero here; the parallel-chord division-by-zero cannot occur.
            let t = point_start / point_end; // plane intersection parameter.

            if !(0.0..=1.0).contains(&t) {
                return None;
            }

            return self.create_intersection(geometry, center, y_start_spatial, direction, t);
        }

        // Same-sign endpoints: the straight chord never crosses the plane, but
        // the true geodesic may dip through it and back within this single step.
        // Reconstruct z(s) as a cubic Hermite from the endpoint heights and
        // z-velocities and look for a hidden crossing.
        //
        // Note: in practice this branch is effectively inert for the black-hole
        // scenes (e.g. Kerr). The adaptive integrator slows for curvature near
        // the equatorial plane, so its steps are too small to straddle a
        // dip-and-back; every real crossing shows up as an endpoint sign change
        // handled by the linear branch above. Instrumenting a full vantage
        // render produced zero hits here (see docs/known-rendering-behaviors.md).
        // The reconstruction is kept because it is correct and cheap and would
        // catch a genuine in-step double crossing if a camera ever produced one.
        // The Hermite slopes are dz/ds with s running over the step's own
        // integration parameter (Step::t). get_z_cartesian returns dz/d(affine
        // lambda); on charts whose Step::t is not the affine parameter this
        // needs rescaling by d(affine)/d(Step::t). KerrBL integrates in Mino
        // time, where d(affine) = Sigma d(Mino), so its slopes carry a factor
        // of Sigma = r^2 + a^2 cos^2(theta); the affine-parametrised charts
        // (Kerr/KS, Schwarzschild) use 1.
        let param_scale = |p: &Point| match p.coordinate_system {
            CoordinateSystem::BoyerLindquist { a } => {
                let (r, theta) = (p[1], p[2]);
                r * r + a * a * theta.cos().powi(2)
            }
            _ => 1.0,
        };
        let d_lambda = y_end.t - y_start.t;
        let p0 = y_start_spatial[2];
        let p1 = y_end_spatial[2];
        let m0 = y_start.p.get_z_cartesian(&y_start.x) * param_scale(&y_start.x) * d_lambda;
        let m1 = y_end.p.get_z_cartesian(&y_end.x) * param_scale(&y_end.x) * d_lambda;

        // p(t) = (2 t^3 - 3 t^2 + 1) p0 + (t^3 - 2 t^2 + t) m0 + (-2 t^3 + 3 t^2) p1 + (t^3 - t^2) m1
        //      = p0 + m0 t + (-2 m0 - m1 - 3 p0 + 3 p1) t^2 + (m0 + m1 + 2 p0 - 2 p1) t^3
        // p'(t) = (6 t^2 - 6 t) p0 + (3 t^2 - 4 t + 1) m0 + (-6 t^2 + 6 t) p1 + (3 t^2 - 2 t) m1
        //       = m0 + (-4 m0 - 2 m1 - 6 p0 + 6 p1) t + (3 m0 + 3 m1 + 6 p0 - 6 p1) t^2

        let c0 = p0;
        let c1 = m0;
        let c2 = -3.0 * p0 - 2.0 * m0 + 3.0 * p1 - m1;
        let c3 = 2.0 * p0 + m0 - 2.0 * p1 + m1;

        let a = 3.0 * c3;
        let b = 2.0 * c2;
        let c = c1;
        let discriminant = b * b - 4.0 * a * c;

        // If the discriminant is negative, there are no real roots, and thus no intersection with the plane.
        // We already checked that the segment does not cross the plane.
        if discriminant < 0.0 {
            return None;
        }
        // Safe way to compute the roots of a quadratic equation. This will not blow up all values
        // if a is very small. See https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
        let q = -0.5 * (b + discriminant.sqrt().copysign(b));
        let int_p1 = q / a;
        let int_p2 = c / q;

        let mut points = vec![0.0, 1.0];
        for tp in [int_p1, int_p2] {
            if tp.is_finite() && tp > 0.0 && tp < 1.0 {
                points.push(tp);
            }
        }
        points.sort_by(f64::total_cmp);

        for window in points.windows(2) {
            let p1 = window[0];
            let p2 = window[1];

            let f = |s: f64| c0 + c1 * s + c2 * s * s + c3 * s * s * s;
            let mut conv = SimpleConvergency {
                eps: 1e-12,
                max_iter: 100,
            };
            let root = find_root_brent(p1, p2, &f, &mut conv);
            let t = match root {
                Ok(t) => t,
                Err(_) => continue,
            };

            let intersection =
                self.create_intersection(geometry, center, y_start_spatial, direction, t);

            if intersection.is_some() {
                return intersection;
            }
        }
        None
    }

    fn color_at_uv(
        &self,
        color_computation_data: &ColorComputationData,
        _geometry: &dyn Geometry,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.texture_mapper.color_at_uv(
            &color_computation_data.uv,
            &color_computation_data.temperature_data,
        )
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

impl SceneObject for Disc {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::kerr::Kerr;
    use crate::rendering::color::Color;
    use crate::rendering::temperature::ConstantTemperatureComputer;
    use crate::rendering::texture::CheckerMapper;
    use std::sync::Arc;

    fn create_disc(inner: f64, outer: f64) -> Disc {
        Disc::new(
            inner,
            outer,
            Arc::new(CheckerMapper::new(
                0.0,
                5.0,
                5.0,
                Color::new(100, 0, 0, 255),
                Color::new(0, 100, 0, 255),
            )),
            Box::new(ConstantTemperatureComputer::new(1000.0)),
        )
    }

    fn crossing_at(x: f64) -> (Step, Step) {
        let start = Point::new_cartesian(0.0, x, 0.0, 0.1);
        let end = Point::new_cartesian(0.0, x, 0.0, -0.1);
        // Straight-line tangent, so the Cubic Hermite reduces to the linear crossing this fixture describes.
        let dir = end.get_spatial_vector_cartesian() - start.get_spatial_vector_cartesian();
        let p = FourVector::new_cartesian(1.0, dir[0], dir[1], dir[2]);
        (
            Step {
                x: start,
                p,
                t: 0.0,
                step: 0,
            },
            Step {
                x: end,
                p,
                t: 1.0,
                step: 1,
            },
        )
    }

    #[test]
    fn test_disc_boundary_uses_metric_radial_coordinate_in_kerr() {
        // Kerr equatorial relation: x^2 + y^2 = r^2 + a^2. A plane crossing
        // at Cartesian radius R corresponds to Boyer-Lindquist
        // r = sqrt(R^2 - a^2), so points slightly outside the inner edge in
        // Cartesian norm can lie INSIDE it in the metric radial coordinate.
        let geometry = Kerr::new(1.0, 0.499, 1e-4);
        let disc = create_disc(0.795, 8.0);

        // Cartesian R = 0.8 -> BL r = sqrt(0.64 - 0.499^2) ~ 0.625 < 0.795:
        // no gas there.
        let (start, end) = crossing_at(0.8);
        assert!(disc.intersects(&start, &end, &geometry).is_none());

        // Cartesian R = 0.95 -> BL r ~ 0.808 > 0.795: inside the disc.
        let (start, end) = crossing_at(0.95);
        assert!(disc.intersects(&start, &end, &geometry).is_some());

        // Outer edge: Cartesian R = 8.01 -> BL r ~ 7.994 < 8.0: still gas.
        let (start, end) = crossing_at(8.01);
        assert!(disc.intersects(&start, &end, &geometry).is_some());
    }

    #[test]
    fn test_disc_boundary_unchanged_in_flat_space() {
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.05, 8.0);

        let (start, end) = crossing_at(3.0);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
        let (start, end) = crossing_at(3.1);
        assert!(disc.intersects(&start, &end, &geometry).is_some());
        let (start, end) = crossing_at(8.1);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
    }

    #[test]
    fn test_disc_uv_normalization_uses_metric_radius() {
        // At the disc's inner edge the normalized UV radius must be ~0
        // (texture inner boundary), evaluated in the metric radial
        // coordinate, not the Cartesian norm.
        let geometry = Kerr::new(1.0, 0.499, 1e-4);
        let disc = create_disc(0.795, 8.0);

        // Cartesian R such that BL r == inner edge: R = sqrt(r^2 + a^2).
        let r_inner_cartesian = (0.795f64 * 0.795 + 0.499 * 0.499).sqrt();
        let (start, end) = crossing_at(r_inner_cartesian + 1e-6);
        let intersection = disc
            .intersects(&start, &end, &geometry)
            .expect("should hit the inner edge");
        // u = 0.5 + 0.5 * r_normalized * cos(phi); at phi = 0 and
        // r_normalized ~ 0 this is ~0.5.
        assert!((intersection.uv.u - 0.5).abs() < 1e-3);
        assert!((intersection.uv.v - 0.5).abs() < 1e-3);
    }

    /// A step with fixed in-plane position (radius `x`, phi = 0), endpoint
    /// heights `z0`/`z1` above the disc plane, and Cartesian z-velocities
    /// `m0`/`m1`. The affine step length is 1, so `get_z_cartesian * d_lambda`
    /// returns `m0`/`m1` directly as the cubic's endpoint slopes. Used to
    /// exercise the same-side cubic branch, where both endpoints share a sign
    /// and the straight chord finds no crossing.
    fn dip_steps(x: f64, z0: f64, z1: f64, m0: f64, m1: f64) -> (Step, Step) {
        (
            Step {
                x: Point::new_cartesian(0.0, x, 0.0, z0),
                p: FourVector::new_cartesian(1.0, 0.0, 0.0, m0),
                t: 0.0,
                step: 0,
            },
            Step {
                x: Point::new_cartesian(0.0, x, 0.0, z1),
                p: FourVector::new_cartesian(1.0, 0.0, 0.0, m1),
                t: 1.0,
                step: 1,
            },
        )
    }

    #[test]
    fn test_cubic_detects_symmetric_same_side_dip() {
        // Both endpoints above the plane (z = 0.1) but the path dips through
        // it and back (m0 < 0, m1 > 0). The straight chord between them never
        // crosses, so only the cubic Hermite catches it. Symmetric slopes make
        // c3 = 0, so this also drives the a -> 0 (parabola) branch of the
        // numerically stable quadratic solve.
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.0, 8.0);
        let (start, end) = dip_steps(5.0, 0.1, 0.1, -1.0, 1.0);
        let hit = disc
            .intersects(&start, &end, &geometry)
            .expect("cubic must catch the same-side dip");
        // The hit must sit on the disc plane, not at the chord's z (0.1);
        // otherwise the radius/UV/emitter are evaluated off the disc.
        assert!(hit.intersection_point.get_spatial_vector_cartesian()[2].abs() < 1e-9);
    }

    #[test]
    fn test_cubic_detects_asymmetric_same_side_dip() {
        // Genuine cubic (c3 != 0): both endpoints positive, asymmetric slopes,
        // still dips below zero and back.
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.0, 8.0);
        let (start, end) = dip_steps(5.0, 0.1, 0.3, -2.0, 1.0);
        assert!(disc.intersects(&start, &end, &geometry).is_some());
    }

    #[test]
    fn test_cubic_same_side_no_dip_misses() {
        // Same-sign endpoints with interior turning points, but the curve
        // never reaches the far side of zero: no crossing, so None. Exercises
        // the "turning points exist but no sub-interval brackets a root" path.
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.0, 8.0);
        let (start, end) = dip_steps(5.0, 0.1, 0.3, 1.0, 1.0);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
    }

    #[test]
    fn test_cubic_same_side_monotone_misses() {
        // Same-sign, strictly monotone: the derivative has no real roots
        // (discriminant < 0), so the early-out returns None without a solve.
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.0, 8.0);
        let (start, end) = dip_steps(5.0, 0.1, 1.0, 1.0, 1.0);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
    }

    #[test]
    fn test_cubic_dip_outside_annulus_misses() {
        // Same dip geometry, but the in-plane radius (1.0) is inside the disc
        // inner edge (3.0), so the crossing lands off the gas: None.
        let geometry = EuclideanSpace::new();
        let disc = create_disc(3.0, 8.0);
        let (start, end) = dip_steps(1.0, 0.1, 0.1, -1.0, 1.0);
        assert!(disc.intersects(&start, &end, &geometry).is_none());
    }
}
