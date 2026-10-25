//! Spherical coordinate frames parameterized by a rotation about the x-axis.

use crate::geometry::spherical_coordinates_helper::rotate_about_x;
use nalgebra::Vector3;

/// Spherical coordinates whose frame is rotated by `rotation` degrees about the x-axis
/// (same convention as `rotate_about_x`: coordinate frame = R(rotation) · physical frame).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct RotatedSpherical {
    pub rotation: f64,
}

pub(in crate::geometry) fn spherical_to_frame(r: f64, theta: f64, phi: f64) -> Vector3<f64> {
    let (st, ct) = theta.sin_cos();
    let (sp, cp) = phi.sin_cos();
    Vector3::new(r * st * cp, r * st * sp, r * ct)
}

pub(in crate::geometry) fn rotate(v: Vector3<f64>, deg: f64) -> Vector3<f64> {
    let (x, y, z) = rotate_about_x(v[0], v[1], v[2], deg);
    Vector3::new(x, y, z)
}
