use crate::geometry::point::Point;
use crate::rendering::color::Color;
use crate::rendering::texture::TextureMap;
use crate::scene_objects::hittable::Hittable;
use nalgebra::RealField;

pub trait SceneObject: Hittable + TextureMap {}

pub struct Objects {
    objects: Vec<Box<dyn SceneObject>>,
}

impl Objects {
    pub fn new() -> Self {
        Self {
            objects: Vec::new(),
        }
    }

    pub fn add_object(&mut self, hittable: Box<dyn SceneObject>) {
        self.objects.push(hittable);
    }

    pub fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<Color> {
        let mut resulting_color = None;
        let mut shortest_distance = f64::MAX;
        let y_start_cartesian = y_start.get_as_cartesian();

        // The step size can be rather large, so it makes sense to sort the objects by their
        // distance to the y_start point.
        for hittable in &self.objects {
            if let Some(uv) = hittable.intersects(y_start, y_end) {
                let intersection = uv.intersection_point.get_as_cartesian();
                let distance = (intersection - y_start_cartesian).norm();
                if distance < shortest_distance {
                    shortest_distance = distance;
                    resulting_color = Some(hittable.color_at_uv(uv.uv));
                }
            }
        }
        resulting_color
    }
}
