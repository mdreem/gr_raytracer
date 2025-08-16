use crate::geometry::point::Point;
use crate::rendering::texture::UVCoordinates;

pub struct Intersection {
    pub uv: UVCoordinates,
    pub intersection_point: Point,
}

pub trait Hittable {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection>;
}
