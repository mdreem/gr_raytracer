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

pub trait TextureMap: Sync {
    fn color_at_uv(&self, uv: UVCoordinates, redshift: f64) -> CIETristimulus;
}

#[derive(Clone)]
pub struct TextureMapper {
    image: DynamicImage,
    color_normalization: CIETristimulusNormalization,
}

impl TextureMapper {
    pub fn new(
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
    fn color_at_uv(&self, uv: UVCoordinates, _redshift: f64) -> CIETristimulus {
        let cie_tristimulus = self.bilinear(&uv);
        cie_tristimulus.normalize(self.color_normalization)
    }
}

#[derive(Clone)]
pub struct BlackBodyMapper {
    temperature: f64,
    color_normalization: CIETristimulusNormalization,
}

impl BlackBodyMapper {
    pub fn new(
        temperature: f64,
        color_normalization: CIETristimulusNormalization,
    ) -> BlackBodyMapper {
        BlackBodyMapper {
            temperature,
            color_normalization,
        }
    }
}

impl TextureMap for BlackBodyMapper {
    fn color_at_uv(&self, _uv: UVCoordinates, redshift: f64) -> CIETristimulus {
        let c = get_cie_xyz_of_black_body_redshifted(self.temperature * redshift);
        let r = c.normalize(self.color_normalization);
        r
    }
}

#[derive(Clone)]
pub struct CheckerMapper {
    width: f64,
    height: f64,
    c1: CIETristimulus,
    c2: CIETristimulus,
    color_normalization: CIETristimulusNormalization,
}

impl CheckerMapper {
    #[allow(dead_code)] // For testing
    pub fn new(
        width: f64,
        height: f64,
        c1: Color,
        c2: Color,
        color_normalization: CIETristimulusNormalization,
    ) -> CheckerMapper {
        CheckerMapper {
            width,
            height,
            c1: srgb_to_xyz(&c1),
            c2: srgb_to_xyz(&c2),
            color_normalization,
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(&self, uv: UVCoordinates, _redshift: f64) -> CIETristimulus {
        let ut = (uv.u * self.width).floor() as usize;
        let vt = (uv.v * self.height).floor() as usize;

        if (ut + vt) % 2 == 0 {
            self.c1.normalize(self.color_normalization)
        } else {
            self.c2.normalize(self.color_normalization)
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
        filename: String,
        color_normalization: CIETristimulusNormalization,
    ) -> Result<TextureMapHandle, RaytracerError> {
        if self.texture_map_map.get(&filename).is_none() {
            let new_mapper = TextureMapper::new(filename.clone(), color_normalization)?;
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
