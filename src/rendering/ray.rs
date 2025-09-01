use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::HasCoordinateSystem;
use crate::geometry::point::Point;
use crate::rendering::integrator::Step;
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::scene::get_position;
use std::io::Write;
use std::ops::Index;

#[derive(Debug)]
pub struct Ray {
    pub position: Point,
    pub momentum: FourVector,
    pub row: i64,
    pub col: i64,
    pub width: i64,
    pub height: i64,
}
impl Ray {
    pub fn new(
        row: i64,
        col: i64,
        width: i64,
        height: i64,
        position: Point,
        momentum: FourVector,
    ) -> Self {
        Self {
            row,
            col,
            width,
            height,
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

    pub fn save(
        &self,
        write: &mut dyn Write,
        geometry: &dyn HasCoordinateSystem,
    ) -> Result<(), RaytracerError> {
        write
            .write_all(b"i,t,tau,x,y,z\n")
            .map_err(RaytracerError::IoError)?;

        for step in &self.steps {
            let position = get_position(&step.y, geometry.coordinate_system());

            write
                .write_all(
                    format!(
                        "{},{},{},{},{},{}\n",
                        step.step, step.t, position[0], position[1], position[2], position[3],
                    )
                    .as_bytes(),
                )
                .map_err(RaytracerError::IoError)?;
        }
        write.flush().map_err(RaytracerError::IoError)
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
