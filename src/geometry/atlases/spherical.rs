//! Rotated spherical charts, their transitions, and pole avoidance.

use super::Atlas;
use crate::geometry::coordinates::chart::Chart;
use crate::geometry::coordinates::spherical::{RotatedSpherical, rotate, spherical_to_frame};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use nalgebra::Vector3;

const POLE_SWITCH_SIN_THETA: f64 = 0.1;

impl Chart for SphericalChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        -1.0
    }
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum SphericalChart {
    Rotated(RotatedSpherical),
}

impl SphericalChart {
    pub fn rotation(&self) -> f64 {
        let SphericalChart::Rotated(RotatedSpherical { rotation }) = self;
        *rotation
    }

    /// The generator of rotations about the physical z-axis (the physical d_phi) at `position`,
    /// in the components of the position's chart. In a rotated chart the chart's own d_phi
    /// rotates about a different axis.
    pub fn physical_axial_vector(position: &Point<Self>) -> FourVector<Self> {
        let rotation = position.chart.rotation();
        if rotation == 0.0 {
            return FourVector::new(0.0, 0.0, 0.0, 1.0, *position);
        }
        let (r, theta, phi) = (position[1], position[2], position[3]);
        let physical = rotate(spherical_to_frame(r, theta, phi), -rotation);
        let d = rotate(Vector3::new(-physical[1], physical[0], 0.0), rotation);

        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        FourVector::new(
            0.0,
            st * cp * d[0] + st * sp * d[1] + ct * d[2],
            (ct * cp * d[0] + ct * sp * d[1] - st * d[2]) / r,
            (-sp * d[0] + cp * d[1]) / (r * st),
            *position,
        )
    }
}

impl Atlas for SphericalChart {
    /// Near the current chart's pole, switch to the chart rotated by 90°, where the point lies
    /// close to the equator. Two perpendicular poles cannot both be near a point
    /// (sin²θ_A + sin²θ_B ≥ 1), so the switch has built-in hysteresis.
    fn best_chart(&self, position: &Point<Self>) -> Self {
        if position[2].sin().abs() >= POLE_SWITCH_SIN_THETA {
            return position.chart;
        }
        SphericalChart::Rotated(RotatedSpherical {
            rotation: (position.chart.rotation() + 90.0).rem_euclid(360.0),
        })
    }

    fn transition_point(&self, position: &Point<Self>, to: Self) -> Point<Self> {
        if position.chart == to {
            return *position;
        }
        let (t, r, theta, phi) = (position[0], position[1], position[2], position[3]);
        let c = rotate(
            spherical_to_frame(r, theta, phi),
            to.rotation() - position.chart.rotation(),
        );
        let theta_to = if r == 0.0 {
            0.0
        } else {
            (c[2] / r).clamp(-1.0, 1.0).acos()
        };
        Point::new(t, r, theta_to, c[1].atan2(c[0]), to)
    }

    fn transition_four_vector(&self, vector: &FourVector<Self>, to: Self) -> FourVector<Self> {
        let at = self.transition_point(&vector.at, to);
        if vector.at.chart == to {
            return *vector;
        }
        let (r, theta, phi) = (vector.at[1], vector.at[2], vector.at[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let (dr, dtheta, dphi) = (vector[1], vector[2], vector[3]);
        let d = rotate(
            Vector3::new(
                st * cp * dr + r * ct * cp * dtheta - r * st * sp * dphi,
                st * sp * dr + r * ct * sp * dtheta + r * st * cp * dphi,
                ct * dr - r * st * dtheta,
            ),
            to.rotation() - vector.at.chart.rotation(),
        );

        let (r, theta, phi) = (at[1], at[2], at[3]);
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        FourVector::new(
            vector[0],
            st * cp * d[0] + st * sp * d[1] + ct * d[2],
            (ct * cp * d[0] + ct * sp * d[1] - st * d[2]) / r,
            (-sp * d[0] + cp * d[1]) / (r * st),
            at,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn rotated(rotation: f64) -> SphericalChart {
        SphericalChart::Rotated(RotatedSpherical { rotation })
    }

    fn spatial_norm(v: &FourVector<SphericalChart>) -> f64 {
        let (r, theta) = (v.at[1], v.at[2]);
        v[1] * v[1] + r * r * v[2] * v[2] + r * r * theta.sin().powi(2) * v[3] * v[3]
    }

    #[test]
    fn rotation_transition_round_trips_points_and_vectors() {
        let p = Point::new(2.0, 5.0, 1.1, 0.7, rotated(0.0));
        let v = FourVector::new(1.0, 0.3, -0.05, 0.02, p);
        let there = p.chart.transition_four_vector(&v, rotated(90.0));
        let back = there.at.chart.transition_four_vector(&there, rotated(0.0));
        assert_abs_diff_eq!(back.at.vector, p.vector, epsilon = 1e-12);
        assert_abs_diff_eq!(back.vector, v.vector, epsilon = 1e-12);
    }

    #[test]
    fn rotation_transition_keeps_the_physical_point_and_the_vector_length() {
        let p = Point::new(0.0, 5.0, 1.1, 0.7, rotated(0.0));
        let q = p.chart.transition_point(&p, rotated(30.0));
        let physical_p = spherical_to_frame(p[1], p[2], p[3]);
        let physical_q = rotate(spherical_to_frame(q[1], q[2], q[3]), -30.0);
        assert_abs_diff_eq!(physical_p, physical_q, epsilon = 1e-12);

        let v = FourVector::new(1.0, 0.3, -0.05, 0.02, p);
        let w = p.chart.transition_four_vector(&v, rotated(30.0));
        assert_abs_diff_eq!(spatial_norm(&v), spatial_norm(&w), epsilon = 1e-12);
        assert_abs_diff_eq!(v[0], w[0]);
    }

    #[test]
    fn best_chart_leaves_a_chart_near_its_pole_only() {
        let away_from_pole = Point::new(0.0, 5.0, 1.0, 0.4, rotated(0.0));
        assert_eq!(
            away_from_pole.chart.best_chart(&away_from_pole),
            rotated(0.0)
        );

        for theta in [0.05, std::f64::consts::PI - 0.05] {
            let near_pole = Point::new(0.0, 5.0, theta, 0.4, rotated(0.0));
            let best = near_pole.chart.best_chart(&near_pole);
            assert_eq!(best, rotated(90.0));
            let switched = near_pole.chart.transition_point(&near_pole, best);
            assert!(switched[2].sin() > 0.99);
            assert_eq!(switched.chart.best_chart(&switched), rotated(90.0));
        }

        let near_pole_of_last = Point::new(0.0, 5.0, 0.05, 0.4, rotated(270.0));
        assert_eq!(
            near_pole_of_last.chart.best_chart(&near_pole_of_last),
            rotated(0.0)
        );
    }

    /// The physical d_phi is the same vector in every chart: the transition of the unrotated
    /// chart's (0, 0, 0, 1).
    #[test]
    fn physical_axial_vector_agrees_across_rotations() {
        let p = Point::new(0.0, 5.0, 1.1, 0.7, rotated(0.0));
        let axial = SphericalChart::physical_axial_vector(&p);
        assert_eq!(axial.vector, FourVector::new(0.0, 0.0, 0.0, 1.0, p).vector);
        for rotation in [30.0, 90.0, 180.0, 270.0] {
            let expected = p.chart.transition_four_vector(&axial, rotated(rotation));
            let actual = SphericalChart::physical_axial_vector(&expected.at);
            assert_abs_diff_eq!(actual.vector, expected.vector, epsilon = 1e-12);
        }
    }

    #[test]
    fn rotating_by_ninety_degrees_moves_the_pole_to_the_equator() {
        let near_pole = Point::new(0.0, 5.0, 1e-3, 0.4, rotated(0.0));
        let q = near_pole.chart.transition_point(&near_pole, rotated(90.0));
        assert_abs_diff_eq!(q[2], std::f64::consts::FRAC_PI_2, epsilon = 1e-2);
    }
}
