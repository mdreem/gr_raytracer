//! Shared texture interfaces and concrete image, blackbody, and checker mappers.

mod blackbody;
mod checker;
mod image;

pub use blackbody::BlackBodyMapper;
pub use checker::CheckerMapper;
pub use image::TextureMapper;

use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::CIETristimulus;
use std::sync::Arc;

#[derive(Debug, thiserror::Error)]
pub enum TextureError {
    #[error("Texture file error: {0}")]
    TextureFileError(String),
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Image decode error: {0}")]
    DecodeError(#[from] ::image::ImageError),
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
    /// Straight XYZ source color, with alpha interpreted by the material:
    /// surface coverage for surfaces, local density mask for volumes.
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError>;
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
