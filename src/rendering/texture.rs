use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::CIETristimulusNormalization::NoNormalization;
use crate::rendering::color::{srgb_to_xyz, CIETristimulus, Color};
use image::{DynamicImage, GenericImageView, ImageReader};
use std::sync::Arc;

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
}

#[derive(Clone)]
pub struct BlackBodyMapper {
    temperature: f64,
}

#[derive(Clone)]
pub struct CheckerMapper {
    width: f64,
    height: f64,
    c1: CIETristimulus,
    c2: CIETristimulus,
}

impl TextureMapper {
    pub fn new(filename: String) -> TextureMapper {
        let image = ImageReader::open(filename)
            .expect("Failed to open image file")
            .decode()
            .expect("Failed to decode image");

        TextureMapper { image }
    }
}

impl TextureMap for TextureMapper {
    fn color_at_uv(&self, uv: UVCoordinates, redshift: f64) -> CIETristimulus {
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
        cie_tristimulus.normalize(NoNormalization)
    }
}

impl TextureMap for BlackBodyMapper {
    fn color_at_uv(&self, _uv: UVCoordinates, redshift: f64) -> CIETristimulus {
        let c = get_cie_xyz_of_black_body_redshifted(self.temperature * redshift);
        let r = c.normalize(NoNormalization);
        r
    }
}

impl CheckerMapper {
    #[allow(dead_code)] // For testing
    pub fn new(width: f64, height: f64, c1: Color, c2: Color) -> CheckerMapper {
        CheckerMapper {
            width,
            height,
            c1: srgb_to_xyz(&c1),
            c2: srgb_to_xyz(&c2),
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(&self, uv: UVCoordinates, _redshift: f64) -> CIETristimulus {
        let ut = (uv.u * self.width).floor() as usize;
        let vt = (uv.v * self.height).floor() as usize;

        if (ut + vt) % 2 == 0 {
            self.c1
        } else {
            self.c2
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

    pub fn get_texture_mapper(&mut self, filename: String) -> TextureMapHandle {
        if self.texture_map_map.get(&filename).is_none() {
            let new_mapper = TextureMapper::new(filename.clone());
            self.texture_map_map
                .insert(filename.clone(), Arc::new(new_mapper));
        }

        self.texture_map_map
            .get(&filename)
            .expect("Failed to get texture mapper")
            .clone()
    }

    pub fn get_blackbody_mapper(&mut self) -> TextureMapHandle {
        let key = "blackbody".to_string();
        if self.texture_map_map.get(&key).is_none() {
            let new_mapper = TextureMapper {
                image: DynamicImage::new_rgb8(1, 1), // Placeholder image
            };
            self.texture_map_map
                .insert(key.clone(), Arc::new(new_mapper));
        }

        self.texture_map_map
            .get(&key)
            .expect("Failed to get blackbody texture mapper")
            .clone()
    }
}
