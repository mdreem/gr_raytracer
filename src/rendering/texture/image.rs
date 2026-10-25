//! Image texture implementation.

use super::{TemperatureData, TextureError, TextureMap, UVCoordinates};
use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::CIETristimulus;
use super::TextureError::DecodeError;
use crate::rendering::light::radiance::{Radiance, RadianceMean};
use image::{DynamicImage, GenericImageView, ImageReader};

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

        // Filter associated light, so RGB hidden in transparent texels cannot
        // bleed into visible edges. Return straight color at this API boundary.
        let mut filtered = RadianceMean::default();
        for (color, weight) in [(c00, w00), (c10, w10), (c01, w01), (c11, w11)] {
            filtered.add(Radiance::from_straight(color), weight);
        }
        filtered
            .mean()
            .unwrap_or(Radiance::TRANSPARENT)
            .to_straight_surface()
    }
}

impl TextureMap for TextureMapper {
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
        self.bilinear(uv)
            .apply_beaming(temperature_data.redshift, self.beaming_exponent)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let color = texture_mapper
            .color_at_uv(
                &uv,
                &TemperatureData {
                    redshift: 1.0,
                    temperature: 0.0,
                },
            )
            .unwrap();
        assert_eq!(color.x, get_red().x);
        assert_eq!(color.y, get_red().y);
        assert_eq!(color.z, get_red().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_bottom_right_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.999, v: 0.999 };
        let color = texture_mapper
            .color_at_uv(
                &uv,
                &TemperatureData {
                    redshift: 1.0,
                    temperature: 0.0,
                },
            )
            .unwrap();
        assert_eq!(color.x, get_red().x);
        assert_eq!(color.y, get_red().y);
        assert_eq!(color.z, get_red().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_bottom_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.0, v: 0.999 };
        let color = texture_mapper
            .color_at_uv(
                &uv,
                &TemperatureData {
                    redshift: 1.0,
                    temperature: 0.0,
                },
            )
            .unwrap();
        assert_eq!(color.x, get_blue().x);
        assert_eq!(color.y, get_blue().y);
        assert_eq!(color.z, get_blue().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_top_right_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.999, v: 0.0 };
        let color = texture_mapper
            .color_at_uv(
                &uv,
                &TemperatureData {
                    redshift: 1.0,
                    temperature: 0.0,
                },
            )
            .unwrap();

        assert_eq!(color.x, get_blue().x);
        assert_eq!(color.y, get_blue().y);
        assert_eq!(color.z, get_blue().z);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn test_texture_mapper_almost_top_left_corner() {
        let texture_mapper = create_texture_mapper();
        let uv = super::UVCoordinates { u: 0.25, v: 0.25 };
        let color = texture_mapper
            .color_at_uv(
                &uv,
                &TemperatureData {
                    redshift: 1.0,
                    temperature: 0.0,
                },
            )
            .unwrap();
        // Premultiplication and conversion back to straight texture color
        // introduce rounding, even when all texels have identical alpha.
        assert_relative_eq!(color.x, (get_red().x + get_blue().x) / 2.0, epsilon = 1e-14);
        assert_relative_eq!(color.y, (get_red().y + get_blue().y) / 2.0, epsilon = 1e-14);
        assert_relative_eq!(color.z, (get_red().z + get_blue().z) / 2.0, epsilon = 1e-14);
        assert_eq!(color.alpha, 128.0 / 255.0);
    }

    #[test]
    fn transparent_texels_do_not_bleed_color_into_filtered_edges() {
        let image = image::RgbaImage::from_fn(2, 2, |x, _| {
            if x == 0 {
                image::Rgba([255, 0, 0, 255])
            } else {
                image::Rgba([0, 0, 255, 0])
            }
        });
        let mapper = TextureMapper {
            beaming_exponent: 0.0,
            image: image::DynamicImage::ImageRgba8(image),
        };
        let edge = mapper.bilinear(&UVCoordinates { u: 0.25, v: 0.25 });
        assert_relative_eq!(edge.as_vector(), get_red().as_vector(), epsilon = 1e-14);
        assert_eq!(edge.alpha, 0.5);
        let transparent = mapper.bilinear(&UVCoordinates { u: 0.5, v: 0.0 });
        assert_eq!(transparent, CIETristimulus::new(0.0, 0.0, 0.0, 0.0));
    }
}
