use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::{srgb_to_xyz, CIETristimulus, CIETristimulusNormalization, Color};
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
}

impl TextureMap for TextureMapper {
    fn color_at_uv(&self, uv: UVCoordinates, _redshift: f64) -> CIETristimulus {
        let (width, height) = self.image.dimensions();
        // If mapping uv to width-1 and height-1 there are distortions horizontally across the
        // center. Maybe because in spherical coordinates u,v < 1.0, i.e. they're not inclusive.
        // Doing it this way seems to stabilize it.
        let pixel = self.image.get_pixel(
            (((width as f64) * uv.u) as u32).min(width - 1),
            (((height as f64) * uv.v) as u32).min(height - 1),
        );
        let mut cie_tristimulus = srgb_to_xyz(&Color::new(pixel[0], pixel[1], pixel[2], pixel[3]));
        cie_tristimulus.alpha = pixel[3] as f64 / 255.0;
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
