use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::{CIETristimulus, Color, srgb_to_xyz};
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::texture::TextureError::DecodeError;
use image::{DynamicImage, GenericImageView, ImageReader};
use log::{error, warn};
use std::sync::Arc;

#[derive(Debug, thiserror::Error)]
pub enum TextureError {
    #[error("Texture file error: {0}")]
    TextureFileError(String),
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Image decode error: {0}")]
    DecodeError(#[from] image::ImageError),
}

pub struct UVCoordinates {
    pub u: f64,
    pub v: f64,
}

pub struct TextureData {
    pub celestial_map: TextureMapHandle,
}

/// Data about the temperature and redshift at a point. Used e.g. for texture mapping with black
/// body radiation.
pub struct TemperatureData {
    pub temperature: f64,
    pub redshift: f64,
}

pub trait TextureMap: Sync {
    fn color_at_uv(&self, uv: &UVCoordinates, temperature_data: &TemperatureData)
    -> CIETristimulus;
}

#[derive(Clone)]
pub struct TextureMapper {
    beaming_exponent: f64,
    image: DynamicImage,
}

impl TextureMapper {
    pub fn new(beaming_exponent: f64, filename: String) -> Result<TextureMapper, RaytracerError> {
        let image = ImageReader::open(filename)
            .map_err(TextureError::IoError)
            .map_err(RaytracerError::TextureError)?
            .decode()
            .map_err(DecodeError)
            .map_err(RaytracerError::TextureError)?;

        Ok(TextureMapper {
            beaming_exponent,
            image,
        })
    }

    /// https://en.wikipedia.org/wiki/Bilinear_interpolation
    fn bilinear(&self, uv: &UVCoordinates) -> CIETristimulus {
        let (width, height) = self.image.dimensions();
        let p_x = (width as f64) * uv.u;
        let p_y = (height as f64) * uv.v;

        // c01 ------- c11
        //  |           |
        //  |           |
        // c00 ------- c10
        let p_x_floor = (p_x.floor() as u32).min(width - 1);
        let p_y_floor = (p_y.floor() as u32).min(height - 1);
        let p_x_ceil = (p_x.ceil() as u32).min(width - 1);
        let p_y_ceil = (p_y.ceil() as u32).min(height - 1);

        let c00 = CIETristimulus::from_rgba(&self.image.get_pixel(p_x_floor, p_y_floor));
        let c01 = CIETristimulus::from_rgba(&self.image.get_pixel(p_x_floor, p_y_ceil));
        let c11 = CIETristimulus::from_rgba(&self.image.get_pixel(p_x_ceil, p_y_ceil));
        let c10 = CIETristimulus::from_rgba(&self.image.get_pixel(p_x_ceil, p_y_floor));

        let dx = p_x - (p_x_floor as f64);
        let dy = p_y - (p_y_floor as f64);

        let w00 = (1.0 - dx) * (1.0 - dy);
        let w01 = (1.0 - dx) * dy;
        let w10 = dx * (1.0 - dy);
        let w11 = dx * dy;

        w00 * c00 + w10 * c10 + w01 * c01 + w11 * c11
    }
}

impl TextureMap for TextureMapper {
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> CIETristimulus {
        self.bilinear(uv)
            .apply_beaming(temperature_data.redshift, self.beaming_exponent)
    }
}

#[derive(Clone)]
pub struct BlackBodyMapper {
    beaming_exponent: f64,
    /// Precomputed blackbody colors indexed by log10(temperature).
    color_profile: Vec<(f64, CIETristimulus)>,
}

const NUM_COLOR_LUT_STEPS: usize = 1000;
const MIN_TEMPERATURE: f64 = 10.0;
const MAX_TEMPERATURE: f64 = 10_000_000.0;

impl BlackBodyMapper {
    pub fn new(beaming_exponent: f64) -> BlackBodyMapper {
        let mut color_profile = Vec::with_capacity(NUM_COLOR_LUT_STEPS);
        let min_log = MIN_TEMPERATURE.log10();
        let max_log = MAX_TEMPERATURE.log10();
        let step = (max_log - min_log) / (NUM_COLOR_LUT_STEPS - 1) as f64;

        for i in 0..NUM_COLOR_LUT_STEPS {
            let log_t = min_log + (i as f64) * step;
            let t = 10.0_f64.powf(log_t);
            let color = get_cie_xyz_of_black_body_redshifted(t, 1.0);
            color_profile.push((log_t, color));
        }

        BlackBodyMapper {
            beaming_exponent,
            color_profile,
        }
    }

    pub fn blackbody_xyz(&self, temperature: f64, redshift: f64) -> CIETristimulus {
        self.sample_blackbody(temperature * redshift)
            .mul_color_part(redshift.powi(-5))
    }

    fn sample_blackbody(&self, temperature: f64) -> CIETristimulus {
        let Some((first_log_t, first_color)) = self.color_profile.first().copied() else {
            error!("BlackBodyMapper color_profile is empty; falling back to black.");
            return CIETristimulus::new(0.0, 0.0, 0.0, 1.0);
        };
        let Some((last_log_t, last_color)) = self.color_profile.last().copied() else {
            error!("BlackBodyMapper color_profile is empty; falling back to black.");
            return CIETristimulus::new(0.0, 0.0, 0.0, 1.0);
        };

        let log_t = temperature.max(MIN_TEMPERATURE).log10();
        if !log_t.is_finite() {
            warn!(
                "BlackBodyMapper received non-finite temperature {}, clamping to LUT minimum.",
                temperature
            );
            return first_color;
        }

        if log_t <= first_log_t {
            return first_color;
        }
        if log_t >= last_log_t {
            return last_color;
        }

        let idx = match self
            .color_profile
            .binary_search_by(|(lt, _)| lt.total_cmp(&log_t))
        {
            Ok(i) => i,
            Err(0) => 0,
            Err(i) => i - 1,
        };

        let (lt0, c0) = self.color_profile[idx];
        let (lt1, c1) = self.color_profile[idx + 1];

        let t = (log_t - lt0) / (lt1 - lt0);

        CIETristimulus::new(
            c0.x + t * (c1.x - c0.x),
            c0.y + t * (c1.y - c0.y),
            c0.z + t * (c1.z - c0.z),
            1.0,
        )
    }
}

impl TextureMap for BlackBodyMapper {
    fn color_at_uv(
        &self,
        _uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> CIETristimulus {
        let redshift = temperature_data.redshift;
        self.blackbody_xyz(temperature_data.temperature, redshift)
            .apply_beaming(redshift, self.beaming_exponent)
    }
}

#[derive(Clone)]
pub struct CheckerMapper {
    beaming_exponent: f64,
    width: f64,
    height: f64,
    c1: CIETristimulus,
    c2: CIETristimulus,
}

impl CheckerMapper {
    #[allow(dead_code)] // For testing
    pub fn new(
        beaming_exponent: f64,
        width: f64,
        height: f64,
        c1: Color,
        c2: Color,
    ) -> CheckerMapper {
        CheckerMapper {
            beaming_exponent,
            width,
            height,
            c1: srgb_to_xyz(&c1),
            c2: srgb_to_xyz(&c2),
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> CIETristimulus {
        let ut = (uv.u * self.width).floor() as usize;
        let vt = (uv.v * self.height).floor() as usize;

        if (ut + vt).is_multiple_of(2) {
            self.c1
                .apply_beaming(temperature_data.redshift, self.beaming_exponent)
        } else {
            self.c2
                .apply_beaming(temperature_data.redshift, self.beaming_exponent)
        }
    }
}

pub type TextureMapHandle = Arc<dyn TextureMap + Send + Sync>;

pub struct TextureMapperFactory {
    texture_map_map: std::collections::HashMap<String, TextureMapHandle>,
}

impl Default for TextureMapperFactory {
    fn default() -> Self {
        Self::new()
    }
}

impl TextureMapperFactory {
    pub fn new() -> Self {
        Self {
            texture_map_map: std::collections::HashMap::new(),
        }
    }

    pub fn get_texture_mapper(
        &mut self,
        beaming_exponent: f64,
        filename: String,
    ) -> Result<TextureMapHandle, RaytracerError> {
        if self.texture_map_map.get(&filename).is_none() {
            let new_mapper = TextureMapper::new(beaming_exponent, filename.clone())?;
            self.texture_map_map
                .insert(filename.clone(), Arc::new(new_mapper));
        }

        self.texture_map_map
            .get(&filename)
            .cloned()
            .ok_or(RaytracerError::TextureError(
                TextureError::TextureFileError(filename.clone()),
            ))
    }
}

#[cfg(test)]
mod tests {
    use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
    use crate::rendering::color::CIETristimulus;
    use crate::rendering::texture::{BlackBodyMapper, TemperatureData, TextureMap, TextureMapper};
    use approx::assert_relative_eq;

    fn get_red() -> CIETristimulus {
        CIETristimulus::from_rgba(&image::Rgba([255, 0, 0, 255]))
    }

    fn get_blue() -> CIETristimulus {
        CIETristimulus::from_rgba(&image::Rgba([0, 0, 255, 255]))
    }

    fn create_texture_mapper() -> TextureMapper {
        let image_buffer = image::ImageBuffer::from_fn(2, 2, |x, y| {
            if (x + y) % 2 == 0 {
                image::Rgba([255, 0, 0, 128])
            } else {
                image::Rgba([0, 0, 255, 128])
            }
        });

        let dynamic_image = image::DynamicImage::ImageRgba8(image_buffer);
        let texture_mapper = TextureMapper {
            beaming_exponent: 3.0,
            image: dynamic_image,
        };
        texture_mapper
    }

    #[test]
    fn test_texture_mapper_top_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.0, v: 0.0 };
        let color = texture_mapper.color_at_uv(
            &uv,
            &TemperatureData {
                redshift: 1.0,
                temperature: 0.0,
            },
        );
        assert_eq!(color.x, get_red().x);
        assert_eq!(color.y, get_red().y);
        assert_eq!(color.z, get_red().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_bottom_right_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.999, v: 0.999 };
        let color = texture_mapper.color_at_uv(
            &uv,
            &TemperatureData {
                redshift: 1.0,
                temperature: 0.0,
            },
        );
        assert_eq!(color.x, get_red().x);
        assert_eq!(color.y, get_red().y);
        assert_eq!(color.z, get_red().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_bottom_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.0, v: 0.999 };
        let color = texture_mapper.color_at_uv(
            &uv,
            &TemperatureData {
                redshift: 1.0,
                temperature: 0.0,
            },
        );
        assert_eq!(color.x, get_blue().x);
        assert_eq!(color.y, get_blue().y);
        assert_eq!(color.z, get_blue().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_top_right_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.999, v: 0.0 };
        let color = texture_mapper.color_at_uv(
            &uv,
            &TemperatureData {
                redshift: 1.0,
                temperature: 0.0,
            },
        );

        assert_eq!(color.x, get_blue().x);
        assert_eq!(color.y, get_blue().y);
        assert_eq!(color.z, get_blue().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_blackbody_lut_matches_direct_integration() {
        let mapper = BlackBodyMapper::new(0.0);

        for &temperature in &[1_000.0, 5_000.0, 10_000.0, 100_000.0] {
            for &redshift in &[0.5, 1.0, 2.0] {
                let lut = mapper.blackbody_xyz(temperature, redshift);
                let direct = get_cie_xyz_of_black_body_redshifted(temperature, redshift);

                assert_relative_eq!(lut.x, direct.x, max_relative = 0.02, epsilon = 1e-14);
                assert_relative_eq!(lut.y, direct.y, max_relative = 0.02, epsilon = 1e-14);
                assert_relative_eq!(lut.z, direct.z, max_relative = 0.02, epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn test_texture_mapper_almost_top_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.25, v: 0.25 };
        let color = texture_mapper.color_at_uv(
            &uv,
            &TemperatureData {
                redshift: 1.0,
                temperature: 0.0,
            },
        );
        assert_eq!(color.x, (get_red().x + get_blue().x) / 2.0);
        assert_eq!(color.y, (get_red().y + get_blue().y) / 2.0);
        assert_eq!(color.z, (get_red().z + get_blue().z) / 2.0);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }
}
