use crate::geometry::coordinates::atlas::{Atlas, Chart};
use crate::geometry::coordinates::stereographic::StereographicChart::{
    WithNorthPole, WithSouthPole,
};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use nalgebra::{Matrix3, Vector3};

const CHART_SWITCH_DISTANCE: f64 = 2.0;

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum StereographicChart {
    WithNorthPole,
    WithSouthPole,
}

impl Chart for StereographicChart {
    // One Chart is the inversion of the other chart, hence the handedness of  the tetrad switches.
    fn tetrad_handedness_fallback(&self) -> f64 {
        match self {
            WithNorthPole => 1.0,
            WithSouthPole => -1.0,
        }
    }
}

fn compute_u_squared<C: Chart>(point: &Point<C>) -> f64 {
    point.vector[1].powi(2) + point.vector[2].powi(2) + point.vector[3].powi(2)
}

impl Atlas for StereographicChart {
    fn best_chart(&self, position: &Point<Self>) -> Self {
        let l_sq = compute_u_squared(position);

        if l_sq <= CHART_SWITCH_DISTANCE * CHART_SWITCH_DISTANCE {
            return position.chart;
        };

        match position.chart {
            WithNorthPole => WithSouthPole,
            WithSouthPole => WithNorthPole,
        }
    }

    fn transition_point(&self, position: &Point<Self>, to: Self) -> Point<Self> {
        if position.chart == to {
            return *position;
        }
        let l_sq = compute_u_squared(position);

        Point::new(
            position.vector[0],
            position.vector[1] / l_sq,
            position.vector[2] / l_sq,
            position.vector[3] / l_sq,
            to,
        )
    }

    fn transition_four_vector(&self, vector: &FourVector<Self>, to: Self) -> FourVector<Self> {
        if vector.at.chart == to {
            return *vector;
        };

        let new_point = self.transition_point(&vector.at, to);
        let j = jacobian(&vector.at);
        let v_spatial = Vector3::new(vector[1], vector[2], vector[3]);
        let new_vector_spatial = j * v_spatial;

        FourVector::new(
            vector[0],
            new_vector_spatial[0],
            new_vector_spatial[1],
            new_vector_spatial[2],
            new_point,
        )
    }
}

fn jacobian(position: &Point<StereographicChart>) -> Matrix3<f64> {
    let (v1, v2, v3) = (position.vector[1], position.vector[2], position.vector[3]);

    let l_sq = v1.powi(2) + v2.powi(2) + v3.powi(2);
    let l_quad = l_sq * l_sq;

    #[rustfmt::skip]
    let data = [
        l_sq - 2.0*v1*v1, - 2.0*v1*v2, - 2.0*v1*v3,
         -2.0*v2*v1, l_sq - 2.0*v2*v2, - 2.0*v2*v3,
        -2.0*v3*v1, -2.0*v3*v2, l_sq - 2.0*v3*v3,
    ];
    Matrix3::from_row_slice(&data) / l_quad
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const NORTH: StereographicChart = StereographicChart::WithNorthPole;
    const SOUTH: StereographicChart = StereographicChart::WithSouthPole;

    fn spatial(v: &FourVector<StereographicChart>) -> Vector3<f64> {
        Vector3::new(v[1], v[2], v[3])
    }

    /// Spatial line element of the unit 3-sphere in stereographic coordinates,
    /// dl^2 = 4 |du|^2 / (1 + |u|^2)^2.
    fn spatial_norm_squared(v: &FourVector<StereographicChart>) -> f64 {
        let u2 = v.at[1].powi(2) + v.at[2].powi(2) + v.at[3].powi(2);
        4.0 * spatial(v).norm_squared() / (1.0 + u2).powi(2)
    }

    #[test]
    fn point_transition_is_an_involution() {
        let p = Point::new(0.7, 1.3, -0.4, 2.1, NORTH);
        let there = NORTH.transition_point(&p, SOUTH);
        assert_eq!(there.chart, SOUTH);
        assert_abs_diff_eq!(there[0], p[0]);
        let back = SOUTH.transition_point(&there, NORTH);
        assert_eq!(back.chart, NORTH);
        assert_abs_diff_eq!(back.vector, p.vector, epsilon = 1e-12);
    }

    #[test]
    fn transition_to_the_same_chart_is_the_identity() {
        let p = Point::new(0.7, 1.3, -0.4, 2.1, NORTH);
        let v = FourVector::new(-1.0, 0.2, 0.5, -0.3, p);
        assert_eq!(NORTH.transition_point(&p, NORTH), p);
        assert_eq!(NORTH.transition_four_vector(&v, NORTH).vector, v.vector);
    }

    #[test]
    fn four_vector_transition_round_trips() {
        let p = Point::new(0.7, 1.3, -0.4, 2.1, NORTH);
        let v = FourVector::new(-1.0, 0.2, 0.5, -0.3, p);
        let there = NORTH.transition_four_vector(&v, SOUTH);
        assert_eq!(there.at.chart, SOUTH);
        assert_abs_diff_eq!(there[0], v[0]);
        let back = SOUTH.transition_four_vector(&there, NORTH);
        assert_abs_diff_eq!(back.at.vector, p.vector, epsilon = 1e-12);
        assert_abs_diff_eq!(back.vector, v.vector, epsilon = 1e-12);
    }

    /// The vector transition must be the derivative of the point transition.
    #[test]
    fn four_vector_transition_matches_finite_differences_of_the_point_map() {
        let p = Point::new(0.7, 1.3, -0.4, 2.1, NORTH);
        let v = FourVector::new(-1.0, 0.2, 0.5, -0.3, p);
        let h = 1e-6;
        let forward = Point::new(
            p[0],
            p[1] + h * v[1],
            p[2] + h * v[2],
            p[3] + h * v[3],
            NORTH,
        );
        let backward = Point::new(
            p[0],
            p[1] - h * v[1],
            p[2] - h * v[2],
            p[3] - h * v[3],
            NORTH,
        );
        let numeric = (NORTH.transition_point(&forward, SOUTH).vector
            - NORTH.transition_point(&backward, SOUTH).vector)
            / (2.0 * h);
        let analytic = NORTH.transition_four_vector(&v, SOUTH);
        assert_abs_diff_eq!(
            Vector3::new(numeric[1], numeric[2], numeric[3]),
            spatial(&analytic),
            epsilon = 1e-8
        );
    }

    /// Both charts describe the same 3-sphere, so the length of a vector must not change.
    #[test]
    fn four_vector_transition_preserves_the_metric() {
        for (p, v) in [
            (Point::new(0.0, 1.3, -0.4, 2.1, NORTH), [0.2, 0.5, -0.3]),
            (Point::new(0.0, 0.1, 0.2, -0.05, NORTH), [1.0, 0.0, 0.0]),
            (Point::new(0.0, -3.0, 4.0, 0.5, SOUTH), [0.0, -0.7, 0.4]),
        ] {
            let v = FourVector::new(-1.0, v[0], v[1], v[2], p);
            let to = if p.chart == NORTH { SOUTH } else { NORTH };
            let w = p.chart.transition_four_vector(&v, to);
            assert_abs_diff_eq!(
                spatial_norm_squared(&w),
                spatial_norm_squared(&v),
                epsilon = 1e-12 * spatial_norm_squared(&v)
            );
        }
    }

    #[test]
    fn jacobian_matches_finite_differences() {
        let p = Point::new(0.0, 1.3, -0.4, 2.1, NORTH);
        let h = 1e-6;
        let j = jacobian(&p);
        for k in 0..3 {
            let mut plus = p;
            let mut minus = p;
            plus.vector[k + 1] += h;
            minus.vector[k + 1] -= h;
            let column = (NORTH.transition_point(&plus, SOUTH).vector
                - NORTH.transition_point(&minus, SOUTH).vector)
                / (2.0 * h);
            for i in 0..3 {
                assert_abs_diff_eq!(j[(i, k)], column[i + 1], epsilon = 1e-8);
            }
        }
    }
    /// The inversion reverses orientation, which is why the two charts carry opposite
    /// handedness fallbacks.
    #[test]
    fn transition_reverses_orientation() {
        for p in [
            Point::new(0.0, 1.3, -0.4, 2.1, NORTH),
            Point::new(0.0, 0.1, 0.2, -0.05, SOUTH),
            Point::new(0.0, -3.0, 4.0, 0.5, NORTH),
        ] {
            assert!(jacobian(&p).determinant() < 0.0, "{p:?}");
        }
        assert_eq!(
            NORTH.tetrad_handedness_fallback(),
            -SOUTH.tetrad_handedness_fallback()
        );
    }

    #[test]
    fn best_chart_switches_only_far_out_and_does_not_flip_back() {
        let near = Point::new(0.0, 0.5, 0.3, -0.2, NORTH);
        assert_eq!(NORTH.best_chart(&near), NORTH);
        for chart in [NORTH, SOUTH] {
            let far = Point::new(0.0, 3.0, -1.0, 0.5, chart);
            let best = chart.best_chart(&far);
            assert_ne!(best, chart);
            let switched = chart.transition_point(&far, best);
            assert_eq!(best.best_chart(&switched), best);
        }
    }
}
