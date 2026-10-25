//! Position and momentum state used by geodesic integration.

use crate::geometry::coordinates::atlas::Chart;
use crate::geometry::point::Point;
use nalgebra::{Const, OVector};

pub type EquationOfMotionState = OVector<f64, Const<8>>;

/// Extract the position in the chart associated with this integration state.
pub fn get_position<C: Chart>(y: &EquationOfMotionState, chart: C) -> Point<C> {
    Point::new(y[0], y[1], y[2], y[3], chart)
}
