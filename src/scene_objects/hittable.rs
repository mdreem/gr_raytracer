use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::color::CIETristimulus;
use crate::rendering::integrator::Step;
use crate::rendering::texture::UVCoordinates;

pub struct Intersection {
    pub uv: UVCoordinates,
    pub intersection_point: Point,
}

pub trait Hittable: Sync {
    fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Intersection>;
    fn color_at_uv(&self, uv: UVCoordinates, redshift: f64) -> CIETristimulus;
    fn energy_of_emitter(&self, geometry: &dyn Geometry, step: &Step) -> f64;
}
