//! Concrete spacetime models and their geodesic solvers.
//!
//! Shared geometry traits live in [`super::geometry`]; coordinate charts live in
//! [`super::coordinates`].

pub mod euclidean;
pub mod euclidean_spherical;
pub mod hypersphere;
pub mod kerr;
pub mod kerr_bl;
pub mod schwarzschild;
