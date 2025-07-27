use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::HasCoordinateSystem;
use crate::rendering::integrator::Step;
use crate::rendering::scene::get_position;
use nalgebra::Vector4;
use std::fs::File;
use std::io::Write;
use std::ops::Index;

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

pub struct IntegratedRay {
    pub steps: Vec<Step>,
}

impl IntegratedRay {
    pub fn new(steps: Vec<Step>) -> Self {
        Self { steps }
    }

    pub fn save(&self, write: &mut dyn Write, geometry: &dyn HasCoordinateSystem) {
        write
            .write_all(b"i,t,tau,x,y,z\n")
            .expect("Unable to write file");

        for step in &self.steps {
            let position = get_position(&step.y, geometry.coordinate_system()).get_as_vector();

            write
                .write_all(
                    format!(
                        "{},{},{},{},{},{}\n",
                        step.step, step.t, position[0], position[1], position[2], position[3],
                    )
                    .as_bytes(),
                )
                .expect("Unable to write file");
        }
    }

    pub fn len(&self) -> usize {
        self.steps.len()
    }
}

impl Index<usize> for IntegratedRay {
    type Output = Step;

    fn index(&self, index: usize) -> &Self::Output {
        &self.steps[index]
    }
}

impl<'a> IntoIterator for &'a IntegratedRay {
    type Item = &'a Step;
    type IntoIter = std::slice::Iter<'a, Step>;

    fn into_iter(self) -> Self::IntoIter {
        self.steps.iter()
    }
}

impl IntegratedRay {
    pub fn iter(&self) -> std::slice::Iter<'_, Step> {
        self.steps.iter()
    }

    pub fn last(&self) -> Option<&Step> {
        self.steps.last()
    }
}
