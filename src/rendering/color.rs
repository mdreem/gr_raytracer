use image::Rgba;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul};

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub alpha: u8,
}

#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Debug)]
pub enum CIETristimulusNormalization {
    NoNormalization,
    Chromaticity,
    EqualLuminance,
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct CIETristimulus {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub alpha: f64,
}

impl CIETristimulus {
    pub fn from_rgba(rgba: &Rgba<u8>) -> CIETristimulus {
        CIETristimulus::from_color(&Color::from_rgba(rgba))
    }

    pub fn from_color(color: &Color) -> CIETristimulus {
        let mut cie_tristimulus = srgb_to_xyz(color);
        cie_tristimulus.alpha = color.alpha as f64 / 255.0;
        cie_tristimulus
    }
    pub fn new(x: f64, y: f64, z: f64, alpha: f64) -> CIETristimulus {
        CIETristimulus { x, y, z, alpha }
    }

    pub fn as_vector(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    pub fn blend(&self, other: &CIETristimulus) -> CIETristimulus {
        // background = self, foreground = other  (“other over self”)
        let ab = self.alpha.clamp(0.0, 1.0);
        let af = other.alpha.clamp(0.0, 1.0);

        let ao = af + ab * (1.0 - af);
        if ao <= 0.0 {
            return CIETristimulus {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                alpha: 0.0,
            };
        }

        let x = (other.x * af + self.x * ab * (1.0 - af)) / ao;
        let y = (other.y * af + self.y * ab * (1.0 - af)) / ao;
        let z = (other.z * af + self.z * ab * (1.0 - af)) / ao;

        CIETristimulus { x, y, z, alpha: ao }
    }

    pub fn normalize(&self, method: CIETristimulusNormalization) -> CIETristimulus {
        match method {
            CIETristimulusNormalization::NoNormalization => *self,
            CIETristimulusNormalization::Chromaticity => {
                let sum = self.x + self.y + self.z;
                if sum == 0.0 {
                    CIETristimulus {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        alpha: self.alpha,
                    }
                } else {
                    CIETristimulus {
                        x: self.x / sum,
                        y: self.y / sum,
                        z: self.z / sum,
                        alpha: self.alpha,
                    }
                }
            }
            CIETristimulusNormalization::EqualLuminance => {
                if self.y == 0.0 {
                    CIETristimulus {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        alpha: self.alpha,
                    }
                } else {
                    let scale = 1.0 / self.y;
                    CIETristimulus {
                        x: self.x * scale,
                        y: 1.0,
                        z: self.z * scale,
                        alpha: self.alpha,
                    }
                }
            }
        }
    }
}

impl Mul<CIETristimulus> for f64 {
    type Output = CIETristimulus;

    fn mul(self, rhs: CIETristimulus) -> Self::Output {
        CIETristimulus {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
            alpha: self * rhs.alpha,
        }
    }
}

impl Add for CIETristimulus {
    type Output = CIETristimulus;

    /// Component-wise addition. Used for bilinear interpolation.
    fn add(self, rhs: Self) -> Self::Output {
        CIETristimulus {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            alpha: self.alpha + rhs.alpha,
        }
    }
}

impl Color {
    pub fn from_rgba(color: &Rgba<u8>) -> Color {
        Color {
            r: color[0],
            g: color[1],
            b: color[2],
            alpha: color[3],
        }
    }
    pub fn new(r: u8, g: u8, b: u8, alpha: u8) -> Color {
        Color { r, g, b, alpha }
    }
}

// https://en.wikipedia.org/wiki/CIE_1931_color_space
fn g(lambda: f64, mu: f64, tau_left: f64, tau_right: f64) -> f64 {
    let tau = if lambda < mu { tau_left } else { tau_right };
    let t = (lambda - mu) * tau;
    (-0.5 * t * t).exp()
}

pub fn x_bar(lambda: f64) -> f64 {
    1.056 * g(lambda, 599.8, 0.0264, 0.0323) + 0.362 * g(lambda, 442.0, 0.0624, 0.0374)
        - 0.065 * g(lambda, 501.1, 0.0490, 0.0382)
}

pub fn y_bar(lambda: f64) -> f64 {
    0.821 * g(lambda, 568.8, 0.0213, 0.0247) + 0.286 * g(lambda, 530.9, 0.0613, 0.0322)
}

pub fn z_bar(lambda: f64) -> f64 {
    1.217 * g(lambda, 437.0, 0.0845, 0.0278) + 0.681 * g(lambda, 459.0, 0.0385, 0.0725)
}

// https://en.wikipedia.org/wiki/SRGB#The_sRGB_transfer_function_.28.22gamma.22.29
fn compand_srgb(linear: f64) -> f64 {
    let sign = if linear < 0.0 { -1.0 } else { 1.0 };
    let a = linear.abs();
    let encoded = if a <= 0.003_130_8 {
        12.92 * a
    } else {
        1.055 * a.powf(1.0 / 2.4) - 0.055
    };
    (sign * encoded).clamp(0.0, 1.0)
}

pub fn xyz_to_linear_srgb(cie_tristimulus: &CIETristimulus) -> Vector3<f64> {
    // 2003 IEC inverse matrix (XYZ -> linear RGB)
    let m = Matrix3::new(
        3.240_625_5,
        -1.537_208_0,
        -0.498_628_6,
        -0.968_930_7,
        1.875_756_1,
        0.041_517_5,
        0.055_710_1,
        -0.204_021_1,
        1.056_995_9,
    );
    m * cie_tristimulus.as_vector()
}

pub fn xyz_to_srgb(cie_tristimulus: &CIETristimulus, exposure: f64) -> Color {
    let mut v_lin = xyz_to_linear_srgb(cie_tristimulus);
    v_lin *= exposure;

    let r = (compand_srgb(v_lin.x.max(0.0)) * 255.0).round() as u8;
    let g = (compand_srgb(v_lin.y.max(0.0)) * 255.0).round() as u8;
    let b = (compand_srgb(v_lin.z.max(0.0)) * 255.0).round() as u8;

    Color {
        r,
        g,
        b,
        alpha: 255,
    }
}

fn inv_compand_srgb(u: f64) -> f64 {
    // u in [0,1] (encoded sRGB) -> linear
    if u <= 0.04045 {
        u / 12.92
    } else {
        ((u + 0.055) / 1.055).powf(2.4)
    }
}

pub fn srgb_to_xyz(color: &Color) -> CIETristimulus {
    let r_s = color.r as f64 / 255.0;
    let g_s = color.g as f64 / 255.0;
    let b_s = color.b as f64 / 255.0;

    let r = inv_compand_srgb(r_s);
    let g = inv_compand_srgb(g_s);
    let b = inv_compand_srgb(b_s);

    let m = Matrix3::new(
        0.412_456_4,
        0.357_576_1,
        0.180_437_5,
        0.212_672_9,
        0.715_152_2,
        0.072_175_0,
        0.019_333_9,
        0.119_192_0,
        0.950_304_1,
    );
    let v_xyz = m * Vector3::new(r, g, b);
    CIETristimulus::new(v_xyz.x, v_xyz.y, v_xyz.z, 1.0)
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_srgb_to_xyz() {
        let color = Color {
            r: 255,
            g: 42,
            b: 10,
            alpha: 255,
        };
        let cie_tristimulus = srgb_to_xyz(&color);
        let color_back = xyz_to_srgb(&cie_tristimulus, 1.0);
        assert_eq!(color, color_back);
    }

    #[test]
    fn cie_blend_retains_background_when_foreground_transparent() {
        let background = CIETristimulus::new(0.2, 0.4, 0.6, 1.0);
        let foreground = CIETristimulus::new(0.8, 0.1, 0.3, 0.0);
        let blended = background.blend(&foreground);

        assert_abs_diff_eq!(blended.x, background.x);
        assert_abs_diff_eq!(blended.y, background.y);
        assert_abs_diff_eq!(blended.z, background.z);
        assert_abs_diff_eq!(blended.alpha, background.alpha);
    }

    #[test]
    fn cie_blend_two_fully_transparent_colors() {
        let background = CIETristimulus::new(0.2, 0.4, 0.6, 0.0);
        let foreground = CIETristimulus::new(0.8, 0.1, 0.3, 0.0);
        let blended = background.blend(&foreground);

        assert_abs_diff_eq!(blended.x, 0.0);
        assert_abs_diff_eq!(blended.y, 0.0);
        assert_abs_diff_eq!(blended.z, 0.0);
        assert_abs_diff_eq!(blended.alpha, 0.0);
    }

    #[test]
    fn cie_blend_mixes_channels() {
        let background = CIETristimulus::new(0.2, 0.4, 0.6, 1.0);
        let foreground = CIETristimulus::new(0.6, 0.4, 0.2, 0.5);
        let blended = background.blend(&foreground);

        assert_abs_diff_eq!(blended.x, 0.4);
        assert_abs_diff_eq!(blended.y, 0.4);
        assert_abs_diff_eq!(blended.z, 0.4);
        assert_abs_diff_eq!(blended.alpha, 1.0);
    }
}
