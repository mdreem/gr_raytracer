use std::f64::consts::E;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub alpha: u8,
}

impl Color {
    pub fn new(r: u8, g: u8, b: u8, alpha: u8) -> Color {
        Color { r, g, b, alpha }
    }

    pub fn get_as_array(&self) -> [u8; 4] {
        [self.r, self.g, self.b, self.alpha]
    }

    pub fn blend(&self, other: &Color) -> Color {
        let alpha_self = self.alpha as f64 / 255.0;
        let alpha_other = other.alpha as f64 / 255.0;

        let alpha_blend = alpha_self + alpha_other * (1.0 - alpha_self);

        Color {
            r: ((self.r as f64 * alpha_self + other.r as f64 * alpha_other * (1.0 - alpha_self))
                / alpha_blend)
                .round() as u8,
            g: ((self.g as f64 * alpha_self + other.g as f64 * alpha_other * (1.0 - alpha_self))
                / alpha_blend)
                .round() as u8,
            b: ((self.b as f64 * alpha_self + other.b as f64 * alpha_other * (1.0 - alpha_self))
                / alpha_blend)
                .round() as u8,
            alpha: (alpha_blend * 255.0).round() as u8,
        }
    }
}

// https://en.wikipedia.org/wiki/CIE_1931_color_space
fn g(x: f64, mu: f64, tau1: f64, tau2: f64) -> f64 {
    let tau = if x < mu { tau1 } else { tau2 };

    E.powf(-tau.powi(2) * (x - mu).powi(2) / 2.0)
}

fn get_cie_xyz(lambda: f64) -> (f64, f64, f64) {
    let x = 1.056 * g(lambda, 599.8, 0.0264, 0.0323) + 0.362 * g(lambda, 442.0, 0.0624, 0.0374)
        - 0.065 * g(lambda, 501.1, 0.0490, 0.0382);
    let y = 0.821 * g(lambda, 568.8, 0.0213, 0.0247) + 0.286 * g(lambda, 530.9, 0.0613, 0.0322);
    let z = 1.217 * g(lambda, 437.0, 0.0845, 0.0278) + 0.681 * g(lambda, 459.0, 0.0385, 0.0725);

    (x, y, z)
}

// https://en.wikipedia.org/wiki/SRGB#The_sRGB_transfer_function_.28.22gamma.22.29
fn xyz_to_srgb(x: f64, y: f64, z: f64) -> Color {
    let r_lin = 3.2406 * x - 1.5372 * y - 0.4986 * z;
    let g_lin = -0.9689 * x + 1.8758 * y + 0.0415 * z;
    let b_lin = 0.0557 * x - 0.2040 * y + 1.0570 * z;

    let r = ((r_lin.clamp(0.0, 1.0)) * 255.0).round() as u8;
    let g = ((g_lin.clamp(0.0, 1.0)) * 255.0).round() as u8;
    let b = ((b_lin.clamp(0.0, 1.0)) * 255.0).round() as u8;

    Color {
        r,
        g,
        b,
        alpha: 255,
    }
}

fn srgb_to_xyz(color: Color) -> (f64, f64, f64) {
    let r = color.r as f64 / 255.0;
    let g = color.g as f64 / 255.0;
    let b = color.b as f64 / 255.0;

    let x = 0.4124564 * r + 0.3575761 * g + 0.1804375 * b;
    let y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
    let z = 0.0193339 * r + 0.1191920 * g + 0.9503041 * b;

    (x, y, z)
}

pub fn wavelength_to_rgb(lambda: f64) -> Color {
    let (x, y, z) = get_cie_xyz(lambda);
    xyz_to_srgb(x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_red_wavelength_to_rgb() {
        let color = wavelength_to_rgb(600.0);
        assert_eq!(
            color,
            Color {
                r: 255,
                g: 42,
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
                r: 44,
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
        let (x, y, z) = srgb_to_xyz(color);
        let color_back = xyz_to_srgb(x, y, z);
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
}
