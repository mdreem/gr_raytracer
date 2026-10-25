use crate::geometry::coordinates::chart::Chart;
use crate::geometry::coordinates::conversion::CartesianConvertible;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::Step;
use std::io::Write;
use std::ops::Index;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray<C: Chart> {
    pub position: Point<C>,
    pub momentum: FourVector<C>,
    pub row: i64,
    pub col: i64,
}
impl<C: Chart> Ray<C> {
    pub fn new(row: i64, col: i64, position: Point<C>, momentum: FourVector<C>) -> Self {
        Self {
            row,
            col,
            position,
            momentum,
        }
    }
}

pub struct IntegratedRay<C: CartesianConvertible> {
    pub steps: Vec<Step<C>>,
}

impl<C: CartesianConvertible> IntegratedRay<C> {
    pub fn new(steps: Vec<Step<C>>) -> Self {
        Self { steps }
    }

    pub fn save(&self, write: &mut dyn Write) -> Result<(), RaytracerError> {
        write
            .write_all(b"i,t,tau,x,y,z\n")
            .map_err(RaytracerError::IoError)?;

        for step in &self.steps {
            let position = step.x.to_cartesian();

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

impl<C: CartesianConvertible> Index<usize> for IntegratedRay<C> {
    type Output = Step<C>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.steps[index]
    }
}

impl<'a, C: CartesianConvertible> IntoIterator for &'a IntegratedRay<C> {
    type Item = &'a Step<C>;
    type IntoIter = std::slice::Iter<'a, Step<C>>;

    fn into_iter(self) -> Self::IntoIter {
        self.steps.iter()
    }
}

impl<C: CartesianConvertible> IntegratedRay<C> {
    pub fn iter(&self) -> std::slice::Iter<'_, Step<C>> {
        self.steps.iter()
    }

    pub fn last(&self) -> Option<&Step<C>> {
        self.steps.last()
    }
}
