//! Chart groupings, exact transitions, and chart selection for each atlas.
//!
//! Coordinate parameters live in [`super::coordinates`]; physical models live in
//! [`super::spacetimes`]. Atlas enums remain the chart identities carried by points
//! and vectors, so a transition preserves the represented event and tangent.

pub mod cartesian;
pub mod kerr;
pub mod spherical;
pub mod stereographic;

use crate::geometry::coordinates::chart::Chart;
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;

/// The chart enum of a geometry is its atlas: it lists every chart, picks the best one for a
/// point and transforms points and vectors between them. Transitions are exact coordinate
/// changes of the same spacetime; the source chart is the one the point or vector carries.
pub trait Atlas: Chart {
    fn best_chart(&self, position: &Point<Self>) -> Self;
    fn transition_point(&self, position: &Point<Self>, to: Self) -> Point<Self>;
    fn transition_four_vector(&self, vector: &FourVector<Self>, to: Self) -> FourVector<Self>;
}
