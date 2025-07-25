use crate::geometry::four_vector::FourVector;
use crate::texture::UVCoordinates;

pub trait Hittable {
    fn intersects(&self, y_start: &FourVector, y_end: &FourVector) -> Option<UVCoordinates>;
}
