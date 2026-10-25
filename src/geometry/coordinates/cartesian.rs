//! Cartesian coordinate chart.

use super::chart::Chart;

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum CartesianChart {
    Cartesian,
}

impl Chart for CartesianChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        -1.0
    }
}
