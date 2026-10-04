use crate::geometry::coordinates::atlas::Chart;

pub struct Cartesian {}

pub struct CartesianAtlas {}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum CartesianChart {
    Cartesian,
}

impl Chart for CartesianChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        1.0
    }
}
