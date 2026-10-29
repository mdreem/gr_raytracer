//! Tone mapping, image conversion, and file output.

use super::Raytracer;
use crate::geometry::geometry::Geometry;
use crate::rendering::error::RaytracerError;
use crate::rendering::light::color::{
    linear_srgb_to_srgb_buffer, xyz_to_linear_srgb, xyz_to_linear_srgb_buffer,
};
use crate::rendering::light::radiance::Radiance;
use image::{ImageBuffer, ImageFormat, Rgb};
use log::info;

/// Clamp a linear channel to what the Radiance RGBE HDR format can store.
/// RGBE shares one 8-bit exponent (bias 128) across the three channels, so a
/// magnitude below 2^-128 underflows the exponent byte; once it wraps, that
/// pixel decodes as a spuriously enormous or infinite value (the disc-edge
/// "fireflies"). Such values are black anyway, so flush them (and negatives)
/// to zero.
fn flush_to_hdr(value: f64) -> f32 {
    const RGBE_MIN_MAGNITUDE: f64 = 2.938_735_877_055_719e-39; // 2^-128
    if value < RGBE_MIN_MAGNITUDE {
        0.0
    } else {
        value as f32
    }
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    /// Render the section to `filename` and return its luminance (CIE Y) per
    /// pixel, row-major, for follow-up checks such as `check_accuracy`.
    pub fn render_section(
        &self,
        from_row: u32,
        from_col: u32,
        to_row: u32,
        to_col: u32,
        filename: String,
    ) -> Result<Vec<f64>, RaytracerError> {
        let luminance: Vec<f64>;
        if filename.ends_with(".hdr") {
            info!("Creating HDR image");
            let raw_cie =
                self.render_section_to_radiance_buffer(from_row, from_col, to_row, to_col)?;
            luminance = raw_cie.iter().map(|c| c.y).collect();
            let buffer: Vec<f32> = raw_cie
                .into_iter()
                .map(|c| xyz_to_linear_srgb(&c.to_xyz_over_black()))
                .flat_map(|c| [flush_to_hdr(c.x), flush_to_hdr(c.y), flush_to_hdr(c.z)])
                .collect();
            let imgbuf_hdr: ImageBuffer<Rgb<f32>, Vec<f32>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf_hdr
                .save_with_format(&filename, ImageFormat::Hdr)
                .map_err(RaytracerError::ImageError)?;
        } else {
            info!("Creating non-HDR image");
            info!(
                "Tone mapping method: {:?}, exposure: {}",
                self.tone_mapping, self.exposure
            );
            let cie_pixels =
                self.render_section_to_radiance_buffer(from_row, from_col, to_row, to_col)?;
            luminance = cie_pixels.iter().map(|c| c.y).collect();
            let cie_pixels = cie_pixels
                .into_iter()
                .map(Radiance::to_xyz_over_black)
                .collect();
            let linear_srgb = xyz_to_linear_srgb_buffer(&cie_pixels);
            let colors = linear_srgb_to_srgb_buffer(&linear_srgb, self.exposure, self.tone_mapping);
            let buffer: Vec<u8> = colors.iter().flat_map(|c| [c.r, c.g, c.b]).collect();
            let imgbuf: ImageBuffer<Rgb<u8>, Vec<u8>> =
                image::ImageBuffer::from_vec(to_col - from_col, to_row - from_row, buffer)
                    .ok_or(RaytracerError::ImageBufferCreation)?;
            imgbuf.save(&filename).map_err(RaytracerError::ImageError)?;
        }

        info!("saved image to {}", filename);
        Ok(luminance)
    }
}

#[cfg(test)]
mod tests {
    use super::flush_to_hdr;

    #[test]
    fn flush_to_hdr_zeroes_rgbe_underflow_but_keeps_representable_values() {
        // Below 2^-128 the shared RGBE exponent byte underflows and wraps,
        // decoding as a spurious huge/inf pixel: those must flush to zero.
        assert_eq!(flush_to_hdr(2.25e-45), 0.0); // the disc-edge "firefly" value
        assert_eq!(flush_to_hdr(1e-121), 0.0);
        assert_eq!(flush_to_hdr(-3.0), 0.0); // negatives clamp too
        // At and above the representable magnitude the value passes through.
        assert_eq!(flush_to_hdr(1e-30), 1e-30_f64 as f32);
        assert_eq!(flush_to_hdr(50_000.0), 50_000.0_f32);
    }
}
