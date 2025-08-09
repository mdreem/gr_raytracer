use crate::geometry::point::Point;
use crate::rendering::texture::UVCoordinates;

pub trait Hittable {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<UVCoordinates>;
}
