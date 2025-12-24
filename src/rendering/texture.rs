use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::{CIETristimulus, CIETristimulusNormalization, Color, srgb_to_xyz};
use crate::rendering::raytracer::RaytracerError;
use crate::rendering::texture::TextureError::DecodeError;
use image::{DynamicImage, GenericImageView, ImageReader};
use std::sync::Arc;

#[derive(Debug)]
pub enum TextureError {
    TextureFileError(String),
    IoError(std::io::Error),
    DecodeError(image::ImageError),
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
    fn color_at_uv(&self, uv: UVCoordinates, temperature_data: TemperatureData) -> CIETristimulus;
}

#[derive(Clone)]
pub struct TextureMapper {
    beaming_exponent: f64,
    image: DynamicImage,
    color_normalization: CIETristimulusNormalization,
}

impl TextureMapper {
    pub fn new(
        beaming_exponent: f64,
        filename: String,
        color_normalization: CIETristimulusNormalization,
    ) -> Result<TextureMapper, RaytracerError> {
        let image = ImageReader::open(filename)
            .map_err(TextureError::IoError)
            .map_err(RaytracerError::TextureError)?
            .decode()
            .map_err(DecodeError)
            .map_err(RaytracerError::TextureError)?;

        Ok(TextureMapper {
            beaming_exponent,
            image,
            color_normalization,
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
    fn color_at_uv(&self, uv: UVCoordinates, temperature_data: TemperatureData) -> CIETristimulus {
        let cie_tristimulus = self.bilinear(&uv);
        let cie_tristimulus_beamed =
            cie_tristimulus.apply_beaming(temperature_data.redshift, self.beaming_exponent);
        cie_tristimulus_beamed.normalize(self.color_normalization)
    }
}

#[derive(Clone)]
pub struct BlackBodyMapper {
    beaming_exponent: f64,
    color_normalization: CIETristimulusNormalization,
}

impl BlackBodyMapper {
    pub fn new(
        beaming_exponent: f64,
        color_normalization: CIETristimulusNormalization,
    ) -> BlackBodyMapper {
        BlackBodyMapper {
            beaming_exponent,
            color_normalization,
        }
    }
}

impl TextureMap for BlackBodyMapper {
    fn color_at_uv(&self, _uv: UVCoordinates, temperature_data: TemperatureData) -> CIETristimulus {
        let c = get_cie_xyz_of_black_body_redshifted(
            temperature_data.temperature * temperature_data.redshift,
        );
        let c_beamed = c.apply_beaming(temperature_data.redshift, self.beaming_exponent);
        let r = c_beamed.normalize(self.color_normalization);
        r
    }
}

#[derive(Clone)]
pub struct CheckerMapper {
    beaming_exponent: f64,
    width: f64,
    height: f64,
    c1: CIETristimulus,
    c2: CIETristimulus,
    color_normalization: CIETristimulusNormalization,
}

impl CheckerMapper {
    #[allow(dead_code)] // For testing
    pub fn new(
        beaming_exponent: f64,
        width: f64,
        height: f64,
        c1: Color,
        c2: Color,
        color_normalization: CIETristimulusNormalization,
    ) -> CheckerMapper {
        CheckerMapper {
            beaming_exponent,
            width,
            height,
            c1: srgb_to_xyz(&c1),
            c2: srgb_to_xyz(&c2),
            color_normalization,
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(&self, uv: UVCoordinates, temperature_data: TemperatureData) -> CIETristimulus {
        let ut = (uv.u * self.width).floor() as usize;
        let vt = (uv.v * self.height).floor() as usize;

        if (ut + vt) % 2 == 0 {
            self.c1
                .apply_beaming(temperature_data.redshift, self.beaming_exponent)
                .normalize(self.color_normalization)
        } else {
            self.c2
                .apply_beaming(temperature_data.redshift, self.beaming_exponent)
                .normalize(self.color_normalization)
        }
    }
}

pub type TextureMapHandle = Arc<dyn TextureMap + Send + Sync>;

pub struct TextureMapperFactory {
    texture_map_map: std::collections::HashMap<String, TextureMapHandle>,
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
        color_normalization: CIETristimulusNormalization,
    ) -> Result<TextureMapHandle, RaytracerError> {
        if self.texture_map_map.get(&filename).is_none() {
            let new_mapper =
                TextureMapper::new(beaming_exponent, filename.clone(), color_normalization)?;
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
    use crate::rendering::color::CIETristimulus;
    use crate::rendering::texture::{TemperatureData, TextureMap, TextureMapper};

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
            color_normalization: super::CIETristimulusNormalization::NoNormalization,
        };
        texture_mapper
    }

    #[test]
    fn test_texture_mapper_top_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.0, v: 0.0 };
        let color = texture_mapper.color_at_uv(
            uv,
            TemperatureData {
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
            uv,
            TemperatureData {
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
            uv,
            TemperatureData {
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
            uv,
            TemperatureData {
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
    fn test_texture_mapper_almost_top_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.25, v: 0.25 };
        let color = texture_mapper.color_at_uv(
            uv,
            TemperatureData {
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
