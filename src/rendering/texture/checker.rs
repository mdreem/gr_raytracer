//! Checker texture implementation.

use super::{TemperatureData, TextureMap, UVCoordinates};
use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::CIETristimulus;
use crate::rendering::light::color::Color;

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
            c1: CIETristimulus::from_color(&c1),
            c2: CIETristimulus::from_color(&c2),
        }
    }
}

impl TextureMap for CheckerMapper {
    fn color_at_uv(
        &self,
        uv: &UVCoordinates,
        temperature_data: &TemperatureData,
    ) -> Result<CIETristimulus, RaytracerError> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checker_texture_preserves_input_alpha() {
        let mapper = super::CheckerMapper::new(
            0.0,
            1.0,
            1.0,
            crate::rendering::light::color::Color::new(255, 0, 0, 64),
            crate::rendering::light::color::Color::new(0, 0, 255, 0),
        );
        let data = TemperatureData {
            temperature: 0.0,
            redshift: 1.0,
        };
        let red = mapper
            .color_at_uv(&UVCoordinates { u: 0.0, v: 0.0 }, &data)
            .unwrap();
        let blue = mapper
            .color_at_uv(&UVCoordinates { u: 1.0, v: 0.0 }, &data)
            .unwrap();
        assert_eq!(red.alpha, 64.0 / 255.0);
        assert_eq!(blue.alpha, 0.0);
    }
}
