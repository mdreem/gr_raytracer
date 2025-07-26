use crate::geometry::four_vector::FourVector;
use crate::rendering::texture::UVCoordinates;

pub trait Hittable {
    fn intersects(&self, y_start: &FourVector, y_end: &FourVector) -> Option<UVCoordinates>;
}
