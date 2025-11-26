use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::CameraError;
use log::debug;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone)]
pub struct Tetrad {
    position: Point,
    pub t: FourVector,
    pub x: FourVector, // vertical wrt. the camera.
    pub y: FourVector, // horizontal wrt. the camera.
    pub z: FourVector, // away from the camera.
}

impl Tetrad {
    pub fn new(
        position: Point,
        t: FourVector,
        x: FourVector,
        y: FourVector,
        z: FourVector,
    ) -> Self {
        Tetrad {
            position,
            t,
            x,
            y,
            z,
        }
    }
}

impl Display for Tetrad {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Tetrad\n")?;
        write!(f, "  position: {:?}\n", self.position)?;
        write!(f, "  t: {:?}\n", self.t)?;
        write!(f, "  x: {:?}\n", self.x)?;
        write!(f, "  y: {:?}\n", self.y)?;
        write!(f, "  z: {:?}\n", self.z)?;
        Ok(())
    }
}

pub struct TetradValidator<'a, G: Geometry> {
    geometry: &'a G,
}

fn is_approx_equal(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

impl<'a, G: Geometry> TetradValidator<'a, G> {
    pub fn new(geometry: &'a G) -> TetradValidator<'a, G> {
        TetradValidator { geometry }
    }

    pub fn validate(&self, tetrad: &Tetrad) -> Result<(), CameraError> {
        let t_t = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.t, &tetrad.t);
        let x_x = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.x, &tetrad.x);
        let y_y = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.y, &tetrad.y);
        let z_z = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.z, &tetrad.z);

        let t_x = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.t, &tetrad.x);
        let t_y = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.t, &tetrad.y);
        let t_z = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.t, &tetrad.z);
        let x_y = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.x, &tetrad.y);
        let x_z = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.x, &tetrad.z);
        let y_z = self
            .geometry
            .inner_product(&tetrad.position, &tetrad.y, &tetrad.z);

        let signature = self.geometry.signature();
        let tol = 1e-5;
        let is_orthonormal = is_approx_equal(t_t, signature[0], tol)
            && is_approx_equal(x_x, signature[1], tol)
            && is_approx_equal(y_y, signature[2], tol)
            && is_approx_equal(z_z, signature[3], tol)
            && is_approx_equal(t_x, 0.0, tol)
            && is_approx_equal(t_y, 0.0, tol)
            && is_approx_equal(t_z, 0.0, tol)
            && is_approx_equal(x_y, 0.0, tol)
            && is_approx_equal(x_z, 0.0, tol)
            && is_approx_equal(y_z, 0.0, tol);

        debug!("position: {:?}", tetrad.position);
        debug!("tetrad: {}", tetrad);
        debug!("inner product checks for tetrad:");
        debug!("  inner product t.t: {:.5}", t_t);
        debug!("  inner product x.x: {:.5}", x_x);
        debug!("  inner product y.y: {:.5}", y_y);
        debug!("  inner product z.z: {:.5}", z_z);
        debug!("");
        debug!("  inner product t.x: {:.5}", t_x);
        debug!("  inner product t.y: {:.5}", t_y);
        debug!("  inner product t.z: {:.5}", t_z);
        debug!("");
        debug!("  inner product x.y: {:.5}", x_y);
        debug!("  inner product x.z: {:.5}", x_z);
        debug!("");
        debug!("  inner product y.z: {:.5}", y_z);
        debug!("");
        debug!("  is_orthonormal: {:?}", is_orthonormal);

        if !is_orthonormal {
            return Err(CameraError::TetradNotOrthonormal);
        }

        Ok(())
    }
}
