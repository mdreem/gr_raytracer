use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub alpha: u8,
}

#[derive(Deserialize, Serialize, Clone, Copy)]
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

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Chromaticity {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Chromaticity {
    pub fn new(x: f64, y: f64, z: f64) -> Chromaticity {
        Chromaticity { x, y, z }
    }

    pub fn as_vector(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }
}

impl Color {
    pub fn new(r: u8, g: u8, b: u8, alpha: u8) -> Color {
        Color { r, g, b, alpha }
    }

    pub fn get_as_array(&self) -> [u8; 4] {
        [self.r, self.g, self.b, self.alpha]
    }

    pub fn blend(&self, other: &Color) -> Color {
        let (r0, g0, b0) = (
            inv_compand_srgb(self.r as f64 / 255.0),
            inv_compand_srgb(self.g as f64 / 255.0),
            inv_compand_srgb(self.b as f64 / 255.0),
        );
        let (r1, g1, b1) = (
            inv_compand_srgb(other.r as f64 / 255.0),
            inv_compand_srgb(other.g as f64 / 255.0),
            inv_compand_srgb(other.b as f64 / 255.0),
        );

        let a0 = self.alpha as f64 / 255.0;
        let a1 = other.alpha as f64 / 255.0;
        let ao = a0 + a1 * (1.0 - a0);
        if ao == 0.0 {
            return Color {
                r: 0,
                g: 0,
                b: 0,
                alpha: 0,
            };
        }

        let r = (r0 * a0 + r1 * a1 * (1.0 - a0)) / ao;
        let g = (g0 * a0 + g1 * a1 * (1.0 - a0)) / ao;
        let b = (b0 * a0 + b1 * a1 * (1.0 - a0)) / ao;

        Color {
            r: (compand_srgb(r).clamp(0.0, 1.0) * 255.0).round() as u8,
            g: (compand_srgb(g).clamp(0.0, 1.0) * 255.0).round() as u8,
            b: (compand_srgb(b).clamp(0.0, 1.0) * 255.0).round() as u8,
            alpha: (ao * 255.0).round() as u8,
        }
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

fn get_cie_xyz(lambda: f64) -> CIETristimulus {
    let x = x_bar(lambda);
    let y = y_bar(lambda);
    let z = z_bar(lambda);

    CIETristimulus::new(x, y, z, 1.0)
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

pub fn wavelength_to_rgb(lambda: f64) -> Color {
    let cie_tristiumlus = get_cie_xyz(lambda);
    xyz_to_srgb(&cie_tristiumlus, 1.0)
}

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_red_wavelength_to_rgb() {
        let color = wavelength_to_rgb(600.0);
        assert_eq!(
            color,
            Color {
                r: 255,
                g: 113,
                b: 0,
                alpha: 255
            }
        );
    }

    #[test]
    fn test_blue_wavelength_to_rgb() {
        let color = wavelength_to_rgb(450.0);
        assert_eq!(
            color,
            Color {
                r: 116,
                g: 0,
                b: 255,
                alpha: 255
            }
        );
    }

    #[test]
    fn test_green_wavelength_to_rgb() {
        let color = wavelength_to_rgb(540.0);
        assert_eq!(
            color,
            Color {
                r: 0,
                g: 255,
                b: 0,
                alpha: 255
            }
        );
    }

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
    fn blend_first_color_stays() {
        let color1 = Color::new(100, 200, 250, 255);
        let color2 = Color::new(0, 0, 0, 255);
        let blended_color = color1.blend(&color2);

        assert_eq!(blended_color.r, 100);
        assert_eq!(blended_color.g, 200);
        assert_eq!(blended_color.b, 250);
        assert_eq!(blended_color.alpha, 255);
    }

    #[test]
    fn blend_two_fully_transparent_colors() {
        let color1 = Color::new(100, 200, 250, 0);
        let color2 = Color::new(0, 0, 0, 0);
        let blended_color = color1.blend(&color2);

        assert_eq!(blended_color.r, 0);
        assert_eq!(blended_color.g, 0);
        assert_eq!(blended_color.b, 0);
        assert_eq!(blended_color.alpha, 0);
    }
}
