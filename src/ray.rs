use crate::geometry::four_vector::FourVector;
use nalgebra::Vector4;

#[derive(Debug)]
pub struct Ray {
    pub position: Vector4<f64>,
    pub momentum: FourVector,
    pub row: i64,
    pub col: i64,
}
impl Ray {
    pub fn new(row: i64, col: i64, position: Vector4<f64>, momentum: FourVector) -> Self {
        Self {
            row,
            col,
            position,
            momentum,
        }
    }
}
