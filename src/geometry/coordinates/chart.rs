//! Chart identity and the capabilities required by a renderable geometry.

use super::conversion::{CartesianConvertible, SphericalConvertible};
use crate::geometry::atlases::Atlas;
use std::fmt::Debug;

pub trait Chart: Copy + PartialEq + Debug + Send + Sync + 'static {
    // TODO: Align all tetrads on one handedness convention (the (phi, theta, r) tetrads are
    // left-handed, the Kerr-Schild Gram-Schmidt one right-handed), then delete this and
    // `spatial_handedness` in camera.rs.
    /// Handedness of the spatial tetrad that geometries build in this chart; the camera
    /// falls back to it where the Cartesian triple product degenerates.
    fn tetrad_handedness_fallback(&self) -> f64;
}

/// Everything the renderer requires of a geometry's chart. Implemented automatically for every
/// chart that meets the bounds; generic code going through `dyn Geometry<Chart = C>` needs
/// `C: GeometryChart`, since the compiler does not infer the bounds of `HasChart::Chart` there.
pub trait GeometryChart: Chart + Atlas + CartesianConvertible + SphericalConvertible {}

impl<T: Chart + Atlas + CartesianConvertible + SphericalConvertible> GeometryChart for T {}
