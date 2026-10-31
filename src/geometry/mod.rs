//! Coordinate charts, atlas transitions, shared geometric tools, and spacetime models.
//!
//! [`coordinates`] defines coordinate parameters and conversion contracts; [`atlases`]
//! groups compatible charts and implements exact transitions and chart selection.
//! [`spacetimes`] defines the physical geometries independently of that organization.

pub mod atlases;
pub mod capture;
pub mod circular_orbit;
pub mod coordinates;
pub mod four_vector;
pub mod geometry;
pub mod gram_schmidt;
pub mod point;
pub mod spacetimes;
pub mod spherical_coordinates_helper;
pub mod tetrad;
