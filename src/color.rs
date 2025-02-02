use std::f64::consts::E;

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
fn xyz_to_srgb(x: f64, y: f64, z: f64) -> (u8, u8, u8) {
    let r_lin = 3.2406 * x - 1.5372 * y - 0.4986 * z;
    let g_lin = -0.9689 * x + 1.8758 * y + 0.0415 * z;
    let b_lin = 0.0557 * x - 0.2040 * y + 1.0570 * z;

    // Apply gamma correction

    let r = ((r_lin.max(0.0).min(1.0)) * 255.0).round() as u8;
    let g = ((g_lin.max(0.0).min(1.0)) * 255.0).round() as u8;
    let b = ((b_lin.max(0.0).min(1.0)) * 255.0).round() as u8;

    (r, g, b)
}

fn wavelength_to_rgb(lambda: f64) -> (u8, u8, u8) {
    let (x, y, z) = get_cie_xyz(lambda);
    xyz_to_srgb(x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_red_wavelength_to_rgb() {
        let (r, g, b) = wavelength_to_rgb(600.0);
        assert_eq!(r, 255);
        assert_eq!(g, 42);
        assert_eq!(b, 0);
    }

    #[test]
    fn test_blue_wavelength_to_rgb() {
        let (r, g, b) = wavelength_to_rgb(450.0);
        assert_eq!(r, 44);
        assert_eq!(g, 0);
        assert_eq!(b, 255);
    }

    #[test]
    fn test_green_wavelength_to_rgb() {
        let (r, g, b) = wavelength_to_rgb(540.0);
        assert_eq!(r, 0);
        assert_eq!(g, 255);
        assert_eq!(b, 0);
    }
}
