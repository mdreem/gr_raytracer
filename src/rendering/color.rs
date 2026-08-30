use crate::rendering::raytracer::RaytracerError;
use clap::ValueEnum;
use image::Rgba;
use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul};
use std::str::FromStr;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub alpha: u8,
}

#[derive(Deserialize, Serialize, Default, Clone, Copy, PartialEq, Debug, ValueEnum)]
pub enum ToneMappingMethod {
    #[default]
    Reinhard,
    GlobalLinear,
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

    /// Applies relativistic beaming effect based on redshift and beaming exponent.
    ///
    /// A nonpositive redshift factor is unphysical (it can reach here from
    /// emitter velocities evaluated outside their validity region) and is
    /// reported as an error: powf would turn it into NaN (fractional
    /// exponent), 1.0 (exponent 0, the default - silently keeping the bogus
    /// color), or infinity (negative exponent).
    pub fn apply_beaming(
        &self,
        redshift: f64,
        beaming_exponent: f64,
    ) -> Result<CIETristimulus, RaytracerError> {
        if redshift <= 0.0 {
            return Err(RaytracerError::UnphysicalRedshift(redshift));
        }
        // z^0 == 1 for every positive z, and a zero exponent is the common
        // case (the physically exact blackbody path uses it), so skip the
        // `pow` call rather than paying for it once per marched sample.
        if beaming_exponent == 0.0 {
            return Ok(*self);
        }
        let beaming_factor = redshift.powf(beaming_exponent);
        Ok(CIETristimulus {
            x: self.x * beaming_factor,
            y: self.y * beaming_factor,
            z: self.z * beaming_factor,
            alpha: self.alpha,
        })
    }

    pub fn mul_color_part(&self, factor: f64) -> CIETristimulus {
        CIETristimulus {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
            alpha: self.alpha,
        }
    }

    pub fn mul_alpha_part(&self, factor: f64) -> CIETristimulus {
        CIETristimulus {
            x: self.x,
            y: self.y,
            z: self.z,
            alpha: self.alpha * factor,
        }
    }

    pub fn mul_all_parts(&self, factor: f64) -> CIETristimulus {
        CIETristimulus {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
            alpha: self.alpha * factor,
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

impl FromStr for Color {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let components = value
            .split(',')
            .map(|component| {
                component.trim().parse::<u8>().map_err(|_| {
                    format!("invalid RGB color '{value}'; expected three values from 0 to 255")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let [r, g, b] = components.as_slice() else {
            return Err(format!(
                "invalid RGB color '{value}'; expected R,G,B (for example 255,0,255)"
            ));
        };
        Ok(Color::new(*r, *g, *b, 255))
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

pub fn xyz_to_linear_srgb_buffer(cie_tristimulus: &Vec<CIETristimulus>) -> Vec<Vector3<f64>> {
    cie_tristimulus.iter().map(xyz_to_linear_srgb).collect()
}

pub fn linear_srgb_to_srgb_buffer(
    linear_rgb: &Vec<Vector3<f64>>,
    exposure: f64,
    tone_mapping_method: ToneMappingMethod,
) -> Vec<Color> {
    // Per-method pre-scale of linear radiance. Reinhard: the user's
    // exposure slides the scene along the L/(1+L) curve. GlobalLinear IS an
    // auto-exposure (normalize the brightest component to white), so the
    // exposure flag deliberately does not apply; multiplying it in would
    // cancel against the normalization anyway.
    let scale = match tone_mapping_method {
        ToneMappingMethod::Reinhard => exposure,
        ToneMappingMethod::GlobalLinear => {
            let max_component = linear_rgb
                .iter()
                .map(|c| c.x.max(c.y).max(c.z))
                .fold(0.0, f64::max);
            if max_component > 0.0 {
                1.0 / max_component
            } else {
                1.0
            }
        }
    };

    linear_rgb
        .iter()
        .map(|c| c * scale)
        .map(|c| match tone_mapping_method {
            ToneMappingMethod::Reinhard => {
                let l_in = 0.2126 * c.x + 0.7152 * c.y + 0.0722 * c.z;
                if l_in > 0.0 {
                    let l_out = l_in / (1.0 + l_in);
                    c * (l_out / l_in)
                } else {
                    c
                }
            }
            ToneMappingMethod::GlobalLinear => c,
        })
        .map(|v_lin| {
            let r = (compand_srgb(v_lin.x.max(0.0)) * 255.0).round() as u8;
            let g = (compand_srgb(v_lin.y.max(0.0)) * 255.0).round() as u8;
            let b = (compand_srgb(v_lin.z.max(0.0)) * 255.0).round() as u8;
            Color {
                r,
                g,
                b,
                alpha: 255,
            }
        })
        .collect()
}

fn inv_compand_srgb(u: f64) -> f64 {
    // u in [0,1] (encoded sRGB) -> linear
    if u <= 0.04045 {
        u / 12.92
    } else {
        ((u + 0.055) / 1.055).powf(2.4)
    }
}

/// Linearized values of all 256 possible 8-bit sRGB channel codes.
///
/// Texture lookups are the hottest consumer of `srgb_to_xyz`: every escaped
/// ray bilinearly samples four texels, and each texel used to cost three
/// `powf(2.4)` calls. The input is an integer code, so the whole transfer
/// function fits in a table and the result is bit-identical to computing it.
static SRGB_TO_LINEAR: std::sync::LazyLock<[f64; 256]> =
    std::sync::LazyLock::new(|| std::array::from_fn(|code| inv_compand_srgb(code as f64 / 255.0)));

pub fn srgb_to_xyz(color: &Color) -> CIETristimulus {
    let r = SRGB_TO_LINEAR[color.r as usize];
    let g = SRGB_TO_LINEAR[color.g as usize];
    let b = SRGB_TO_LINEAR[color.b as usize];

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

    /// Exposure scales linear radiance before the tone map: Reinhard
    /// responds (dim pixels brighten ~linearly), while GlobalLinear's
    /// auto-normalization cancels it exactly (it is an auto-exposure;
    /// documented on the --exposure flag).
    #[test]
    fn test_exposure_scales_reinhard_but_cancels_in_global_linear() {
        let buffer = vec![Vector3::new(0.02, 0.02, 0.02), Vector3::new(0.5, 0.5, 0.5)];

        let reinhard_1x = linear_srgb_to_srgb_buffer(&buffer, 1.0, ToneMappingMethod::Reinhard);
        let reinhard_4x = linear_srgb_to_srgb_buffer(&buffer, 4.0, ToneMappingMethod::Reinhard);
        assert!(reinhard_4x[0].r > reinhard_1x[0].r);

        let linear_1x = linear_srgb_to_srgb_buffer(&buffer, 1.0, ToneMappingMethod::GlobalLinear);
        let linear_4x = linear_srgb_to_srgb_buffer(&buffer, 4.0, ToneMappingMethod::GlobalLinear);
        assert_eq!(linear_1x, linear_4x);
    }

    /// An unphysical nonpositive redshift must error for every exponent
    /// regime: fractional (would be NaN via powf), zero (would be 1.0,
    /// silently keeping the bogus color), and negative (would be infinity).
    #[test]
    fn test_apply_beaming_rejects_nonpositive_redshift() {
        let c = CIETristimulus::new(1.0, 1.0, 1.0, 1.0);
        for exponent in [3.5, 0.0, -2.0] {
            for redshift in [-0.5, 0.0] {
                assert!(
                    matches!(
                        c.apply_beaming(redshift, exponent),
                        Err(RaytracerError::UnphysicalRedshift(_))
                    ),
                    "exp {} z {}",
                    exponent,
                    redshift
                );
            }
        }
        // Physical redshift still beams normally.
        assert!(c.apply_beaming(0.5, 3.0).unwrap().x > 0.0);
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
    fn parses_rgb_color_with_optional_whitespace() {
        assert_eq!(
            " 12, 34 ,56 ".parse::<Color>().unwrap(),
            Color::new(12, 34, 56, 255)
        );
    }

    #[test]
    fn rejects_malformed_rgb_colors() {
        for value in ["12,34", "12,34,56,78", "12,blue,56", "12,34,256"] {
            assert!(
                value.parse::<Color>().is_err(),
                "{value} should be rejected"
            );
        }
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
