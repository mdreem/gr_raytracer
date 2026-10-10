use crate::geometry::coordinates::atlas::{Atlas, Chart};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum CartesianChart {
    Cartesian,
}

impl Chart for CartesianChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        -1.0
    }
}

impl Atlas for CartesianChart {
    fn best_chart(&self, position: &Point<Self>) -> Self {
        position.chart
    }

    fn transition_point(&self, position: &Point<Self>, _to: Self) -> Point<Self> {
        *position
    }

    fn transition_four_vector(&self, vector: &FourVector<Self>, _to: Self) -> FourVector<Self> {
        *vector
    }
}
