//! Kerr-Schild coordinate parameters and radial/angular helpers.

use nalgebra::Vector4;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct KerrSchildData {
    pub a: f64,
    pub r_s: f64,
}

impl KerrSchildData {
    /// BL r of a point given in Kerr-Schild coordinates (t, x, y, z).
    pub fn bl_radius(&self, v: &Vector4<f64>) -> f64 {
        let (x, y, z) = (v[1], v[2], v[3]);
        let rho2 = x * x + y * y + z * z;
        let a2 = self.a * self.a;
        (0.5 * (rho2 - a2 + ((rho2 - a2).powi(2) + 4.0 * a2 * z * z).sqrt())).sqrt()
    }

    /// sin theta_BL of a point given in Kerr-Schild coordinates (t, x, y, z): cos theta = z / r.
    pub fn sin_theta(&self, v: &Vector4<f64>) -> f64 {
        let r = self.bl_radius(v);
        if r == 0.0 {
            return 0.0;
        }
        let cos_theta = (v[3] / r).clamp(-1.0, 1.0);
        (1.0 - cos_theta * cos_theta).sqrt()
    }
}
