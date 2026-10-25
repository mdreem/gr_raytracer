//! Identity transitions for the single Cartesian chart.

use super::Atlas;
use crate::geometry::coordinates::cartesian::CartesianChart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;

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
