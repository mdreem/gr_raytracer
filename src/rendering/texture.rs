use crate::rendering::color::Color;
use image::{DynamicImage, GenericImageView, ImageReader};

pub struct UVCoordinates {
    pub u: f64,
    pub v: f64,
}

pub struct TextureData<T: TextureMap> {
    pub celestial_map: T,
    pub center_disk_map: T,
    pub center_sphere_map: T,
}

pub trait TextureMap: Sync {
    fn color_at_uv(&self, uv: UVCoordinates) -> Color;
}

#[derive(Clone)]
pub struct TextureMapper {
    image: DynamicImage,
}

#[derive(Clone)]
pub struct CheckerMapper {
    width: f64,
    height: f64,
    c1: Color,
    c2: Color,
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
    fn color_at_uv(&self, uv: UVCoordinates) -> Color {
        let (width, height) = self.image.dimensions();
        let pixel = self.image.get_pixel(
            (((width - 1) as f64) * uv.u) as u32,
            (((height - 1) as f64) * uv.v) as u32,
        );
        Color::new(pixel[0], pixel[1], pixel[2])
    }
}

impl CheckerMapper {
    pub fn new(width: f64, height: f64, c1: Color, c2: Color) -> CheckerMapper {
        CheckerMapper {
            width,
            height,
            c1,
            c2,
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(&self, uv: UVCoordinates) -> Color {
        let ut = (uv.u * self.width).floor() as usize;
        let vt = (uv.v * self.height).floor() as usize;

        if (ut + vt) % 2 == 0 {
            self.c1
        } else {
            self.c2
        }
    }
}
