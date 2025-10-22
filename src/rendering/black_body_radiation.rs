use crate::rendering::color::{CIETristimulus, x_bar, y_bar, z_bar};

const PLANCK_CONSTANT: f64 = 6.626_070_15e-34;
const SPEED_OF_LIGHT: f64 = 299_792_458.0;
const BOLTZMANN_CONSTANT: f64 = 1.380_649e-23;

const MIN_WAVELENGTH: f64 = 380.0; // nm
const MAX_WAVELENGTH: f64 = 830.0; // nm
const NM_TO_M: f64 = 1e-9;

// Planck law is written using meters.
fn planck_spectral_radiance(lambda: f64, temperature: f64) -> f64 {
    let a = 2.0 * PLANCK_CONSTANT * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let b = PLANCK_CONSTANT * SPEED_OF_LIGHT / (lambda * BOLTZMANN_CONSTANT * temperature);
    a / (lambda.powi(5) * (b.exp() - 1.0))
}

fn integrate_blackbody_xyz(temperature: f64) -> CIETristimulus {
    let interval = (MAX_WAVELENGTH - MIN_WAVELENGTH) * NM_TO_M; // in meters
    let step_size = 1.0 * NM_TO_M; // in meters
    let num_steps = ((interval / step_size).floor()) as usize;

    let mut x_accum = 0.0;
    let mut y_accum = 0.0;
    let mut z_accum = 0.0;

    for i in 0..(num_steps as usize) {
        let lambda = MIN_WAVELENGTH * NM_TO_M + (i as f64 + 0.5) * step_size;
        let radiance = planck_spectral_radiance(lambda, temperature);
        x_accum += radiance * x_bar(lambda / NM_TO_M) * step_size;
        y_accum += radiance * y_bar(lambda / NM_TO_M) * step_size;
        z_accum += radiance * z_bar(lambda / NM_TO_M) * step_size;
    }
    CIETristimulus::new(x_accum, y_accum, z_accum, 1.0)
}

pub fn get_cie_xyz_of_black_body_redshifted(temperature: f64) -> CIETristimulus {
    integrate_blackbody_xyz(temperature)
}

#[cfg(test)]
mod tests {
    use super::integrate_blackbody_xyz;
    use crate::rendering::color::{Color, xyz_to_linear_srgb, xyz_to_srgb};
    use image::{ImageFormat, Rgb};

    fn get_srgb_of_black_body(temperature: f64) -> Color {
        get_srgb_of_black_body_redshifted(temperature, 1.0)
    }

    fn get_srgb_of_black_body_redshifted(temperature: f64, redshift: f64) -> Color {
        let cie_tristimulus = integrate_blackbody_xyz(temperature * redshift);
        let exposure = 1.0 / (cie_tristimulus.x + cie_tristimulus.y + cie_tristimulus.z);
        xyz_to_srgb(&cie_tristimulus, exposure)
    }

    #[test]
    fn test_black_body_radiation_red() {
        let color = get_srgb_of_black_body(1000.0);
        assert_eq!(color, Color::new(255, 60, 0, 255));
    }

    #[test]
    fn test_black_body_radiation_blue() {
        let color = get_srgb_of_black_body(10000.0);
        assert_eq!(color, Color::new(137, 146, 172, 255));
    }

    #[test]
    #[ignore]
    fn write_black_body_spectrum_to_file_hdr() {
        let mut imgbuf = image::ImageBuffer::new(100, 100);
        let max_temp_steps = 100;
        let max_redshift_steps = 100;

        for x in 0..max_temp_steps {
            for y in 0..max_redshift_steps {
                let redshift = 0.5 + (y as f64) * 2.0 / (max_temp_steps as f64 - 1.0);
                let temperature = 1000.0 + (x as f64) * 10000.0 / (max_temp_steps as f64 - 1.0);
                let color_xyz = integrate_blackbody_xyz(temperature * redshift);
                let color = xyz_to_linear_srgb(&color_xyz);
                imgbuf.put_pixel(x, y, Rgb([color.x as f32, color.y as f32, color.z as f32]));
            }
        }

        let filename = "black_body_spectrum.hdr";
        imgbuf.save_with_format(filename, ImageFormat::Hdr).unwrap();
        println!("saved image to {}", filename);
    }

    #[test]
    #[ignore]
    fn write_black_body_spectrum_to_file() {
        let mut imgbuf = image::ImageBuffer::new(1000, 1000);
        let max_temp_steps = 100;
        let max_redshift_steps = 100;
        for x in 0..max_temp_steps {
            for y in 0..max_redshift_steps {
                let redshift = 0.5 + (y as f64) * 2.0 / (max_temp_steps as f64 - 1.0);
                let temperature = 1000.0 + (x as f64) * 10000.0 / (max_temp_steps as f64 - 1.0);
                let color = get_srgb_of_black_body_redshifted(temperature, redshift);
                imgbuf.put_pixel(x, y, image::Rgba([color.r, color.g, color.b, color.alpha]));
            }
            println!("line: {}", x)
        }

        let filename = "black_body_spectrum.png";
        imgbuf.save(filename).unwrap();
        println!("saved image to {}", filename);
    }
}
