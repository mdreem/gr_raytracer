use crate::geometry::geometry::Geometry;
use crate::rendering::color::CIETristimulus;
use crate::rendering::redshift::RedshiftComputer;
use crate::rendering::scene::{get_position, EquationOfMotionState};
use crate::scene_objects::hittable::Hittable;

pub trait SceneObject: Hittable {}

pub struct Objects<'a, G: Geometry> {
    geometry: &'a G,
    objects: Vec<Box<dyn SceneObject>>,
}

impl<'a, G: Geometry> Objects<'a, G> {
    pub fn new(geometry: &'a G) -> Self {
        Self {
            geometry,
            objects: Vec::new(),
        }
    }

    pub fn add_object(&mut self, hittable: Box<dyn SceneObject>) {
        self.objects.push(hittable);
    }

    pub fn intersects(
        &self,
        y_start: &EquationOfMotionState,
        y_end: &EquationOfMotionState,
        observer_energy: f64,
        redshift_computer: &RedshiftComputer<'a, G>,
    ) -> Option<CIETristimulus> {
        let mut resulting_color = None;
        let mut shortest_distance = f64::MAX;

        let y_start_point = get_position(&y_start, self.geometry.coordinate_system());
        let y_end_point = get_position(&y_end, self.geometry.coordinate_system());
        let y_start_cartesian = y_start_point.get_as_cartesian();

        // The step size can be rather large, so it makes sense to sort the objects by their
        // distance to the y_start point.
        for hittable in &self.objects {
            if let Some(intersection_data) = hittable.intersects(&y_start_point, &y_end_point) {
                let intersection_point = intersection_data.intersection_point.get_as_cartesian();
                let distance = (intersection_point - y_start_cartesian).norm();
                if distance < shortest_distance {
                    shortest_distance = distance;
                    let redshift = redshift_computer.compute_redshift(y_start, observer_energy);
                    resulting_color = Some(hittable.color_at_uv(intersection_data.uv, redshift));
                }
            }
        }
        resulting_color
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::point::Point;
    use crate::rendering::color::CIETristimulusNormalization::NoNormalization;
    use crate::rendering::color::Color;
    use crate::rendering::texture::CheckerMapper;
    use crate::scene_objects::sphere::Sphere;
    use approx::assert_abs_diff_eq;
    use std::sync::Arc;

    fn create_sphere_at(x: f64, y: f64, z: f64, radius: f64, color: u8) -> Box<dyn SceneObject> {
        Box::new(Sphere::new(
            radius,
            Arc::new(CheckerMapper::new(
                5.0,
                5.0,
                Color::new(color, color, color, 255),
                Color::new(color, color, color, 255),
                NoNormalization,
            )),
            Point::new_cartesian(0.0, x, y, z),
        ))
    }
    #[test]
    #[ignore] // TODO: if step goes through sphere fully, not intersection seen. Need to fix.
    fn test_add_and_intersect_two_sphere() {
        let geometry = EuclideanSpace::new();
        let mut objects = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 1.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects.add_object(farther_sphere);
        objects.add_object(closer_sphere);

        let y_start =
            EquationOfMotionState::from_column_slice(&[0.0, 0.0, 0.0, -3.0, 1.0, 0.0, 0.0, 0.0]);
        let y_end =
            EquationOfMotionState::from_column_slice(&[0.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0, 0.0]);

        let result = objects.intersects(&y_start, &y_end, 1.0, &RedshiftComputer::new(&geometry));
        assert!(result.is_some());
    }

    #[test]
    fn test_add_and_intersect_spheres_inside_each_other() {
        let y_start =
            EquationOfMotionState::from_column_slice(&[0.0, 0.0, 0.0, -3.0, 1.0, 0.0, 0.0, 0.0]);
        let y_end =
            EquationOfMotionState::from_column_slice(&[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);

        let geometry = EuclideanSpace::new();
        let mut objects_setup_1 = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);

        objects_setup_1.add_object(farther_sphere);
        objects_setup_1.add_object(closer_sphere);

        let result_1 =
            objects_setup_1.intersects(&y_start, &y_end, 1.0, &RedshiftComputer::new(&geometry));
        assert!(result_1.is_some());
        assert_abs_diff_eq!(result_1.unwrap().x, 0.121, epsilon = 1e-2);

        let mut objects_setup_2 = Objects::new(&geometry);
        let closer_sphere = create_sphere_at(0.0, 0.0, 0.0, 2.0, 100);
        let farther_sphere = create_sphere_at(0.0, 0.0, 1.0, 1.0, 200);
        objects_setup_2.add_object(closer_sphere);
        objects_setup_2.add_object(farther_sphere);

        let result_2 =
            objects_setup_2.intersects(&y_start, &y_end, 1.0, &RedshiftComputer::new(&geometry));
        assert!(result_2.is_some());
        assert_abs_diff_eq!(result_2.unwrap().x, 0.121, epsilon = 1e-2);
    }
}
