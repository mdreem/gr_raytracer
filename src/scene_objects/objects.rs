use crate::geometry::four_vector::FourVector;
use crate::rendering::color::Color;
use crate::rendering::texture::TextureMap;
use crate::scene_objects::hittable::Hittable;

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

    pub fn intersects(&self, y_start: &FourVector, y_end: &FourVector) -> Option<Color> {
        for hittable in &self.objects {
            if let Some(uv) = hittable.intersects(y_start, y_end) {
                return Some(hittable.color_at_uv(uv));
            }
        }
        None
    }
}
