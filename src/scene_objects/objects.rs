use crate::geometry::point::Point;
use crate::rendering::color::{CIETristimulus, Color};
use crate::rendering::texture::TextureMap;
use crate::scene_objects::hittable::Hittable;
use std::sync::Arc;

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

    pub fn intersects(&self, y_start: &Point, y_end: &Point) -> Option<CIETristimulus> {
        let mut resulting_color = None;
        let mut shortest_distance = f64::MAX;
        let y_start_cartesian = y_start.get_as_cartesian();

        // The step size can be rather large, so it makes sense to sort the objects by their
        // distance to the y_start point.
        for hittable in &self.objects {
            if let Some(intersection_data) = hittable.intersects(y_start, y_end) {
                let intersection_point = intersection_data.intersection_point.get_as_cartesian();
                let distance = (intersection_point - y_start_cartesian).norm();
                if distance < shortest_distance {
                    shortest_distance = distance;
                    resulting_color = Some(hittable.color_at_uv(intersection_data.uv));
                }
            }
        }
        resulting_color
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::point::Point;
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects::hittable::Hittable;
    use crate::scene_objects::sphere::Sphere;
    use std::sync::Arc;

    fn create_sphere_at(x: f64, y: f64, z: f64, radius: f64, color: u8) -> Box<dyn SceneObject> {
        Box::new(Sphere::new(
            radius,
            Arc::new(CheckerMapper::new(
                5.0,
                5.0,
                Color::new(color, color, color, 255),
                Color::new(color, color, color, 255),
            )),
            Point::new_cartesian(0.0, x, y, z),
        ))
    }
    #[test]
    #[ignore] // TODO: if step goes through sphere fully, not intersection seen. Need to fix.
    fn test_add_and_intersect_two_sphere() {
        let mut objects = Objects::new();
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 1.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects.add_object(farther_sphere);
        objects.add_object(closer_sphere);

        let y_start = Point::new_cartesian(0.0, 0.0, 0.0, -3.0);
        let y_end = Point::new_cartesian(0.0, 0.0, 0.0, 3.0);

        let result = objects.intersects(&y_start, &y_end);
        assert!(objects.intersects(&y_start, &y_end).is_some());
    }

    #[test]
    fn test_add_and_intersect_spheres_inside_each_other() {
        let y_start = Point::new_cartesian(0.0, 0.0, 0.0, -3.0);
        let y_end = Point::new_cartesian(0.0, 0.0, 0.0, 0.0);

        let mut objects_setup_1 = Objects::new();
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects_setup_1.add_object(farther_sphere);
        objects_setup_1.add_object(closer_sphere);

        let result_1 = objects_setup_1.intersects(&y_start, &y_end);
        assert!(result_1.is_some());
        assert_eq!(result_1.unwrap().r, 100);

        let mut objects_setup_2 = Objects::new();
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);
        objects_setup_2.add_object(closer_sphere);
        objects_setup_2.add_object(farther_sphere);

        let result_2 = objects_setup_2.intersects(&y_start, &y_end);
        assert!(result_2.is_some());
        assert_eq!(result_2.unwrap().r, 100);
    }
}
