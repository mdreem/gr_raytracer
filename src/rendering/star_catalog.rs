//! Loader for the Gaia DR3 star-catalogue subset downloaded by
//! `scripts/gaia/download.py` (a Zstandard-compressed Parquet file).
//!
//! This is the renderer-side entry point of the point-source star pipeline
//! (see `docs/plan-12-star-catalog.md`). It reads the catalogue columns and
//! precomputes, per star, a unit direction on the celestial sphere. Working in
//! cartesian unit vectors keeps the later gather free of the (theta, phi)
//! coordinate singularities (pole and phi-seam): every membership / solid-angle
//! test is then a plain dot / cross product.

use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::CIETristimulus;
use arrow::array::{Float32Array, Float64Array, Int64Array};
use nalgebra::Vector3;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::file::reader::ChunkReader;
use std::fs::File;
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum StarCatalogError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Parquet error: {0}")]
    ParquetError(#[from] parquet::errors::ParquetError),
    #[error("Arrow error: {0}")]
    ArrowError(#[from] arrow::error::ArrowError),
    #[error("Missing expected column {0:?}")]
    MissingColumn(&'static str),
    #[error("Column {0:?} has unexpected type {1:?}")]
    UnexpectedColumnType(&'static str, arrow::datatypes::DataType),
}

/// A single catalogue star.
#[derive(Debug, Clone)]
pub struct Star {
    pub source_id: i64,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub g_mag: f64,
    pub bp_mag: f64,
    pub rp_mag: f64,
    pub bp_rp: f64,
    /// Unit direction on the celestial sphere (ICRS), precomputed from ra/dec.
    pub direction: Vector3<f64>,
    /// Colour temperature (K) derived from `bp_rp` (see
    /// `colour_temperature_from_bp_rp`). Kept so a later redshifted re-tint
    /// (`T_obs = z * T`) can recompute the chroma per tube.
    pub temperature: f64,
    /// Per-star XYZ radiance the gather sums: the G-band flux (Pogson) times
    /// the star's Y-normalised blackbody chromaticity, so the vector's
    /// luminance (Y) equals the observed flux while X and Z carry its hue.
    /// Precomputed once because it depends only on `g_mag` and `bp_rp`.
    pub emission_xyz: CIETristimulus,
}

impl Star {
    /// Relative linear flux from the G magnitude (Pogson's law), normalized so
    /// magnitude 0 maps to flux 1. Absolute calibration (magnitude 0 -> chosen
    /// linear luminance at unit exposure) is a later concern of the gather.
    pub fn relative_flux(&self) -> f64 {
        pogson_flux(self.g_mag)
    }
}

/// Relative linear flux from a magnitude (Pogson's law): magnitude 0 -> 1,
/// each 5 magnitudes a factor of 100.
fn pogson_flux(magnitude: f64) -> f64 {
    10f64.powf(-0.4 * magnitude)
}

/// Colour temperature (K) from the Gaia BP-RP colour index.
///
/// Uses a Ballesteros-style reciprocal `T = a / (bp_rp + b) + c`. The
/// reciprocal form (rather than a plain polynomial in `bp_rp`) follows from
/// Wien's law, where the colour index runs roughly as `1/T`; a polynomial fit
/// misbehaves at the blue and red ends. The three constants are fit to
/// reference `(BP-RP, Teff)` points spanning the sequence: the hot end
/// (`bp_rp = 0 -> ~9500 K`), the Sun (`0.82 -> 5772 K`), and the cool end
/// (`3.0 -> ~3300 K`), taken from the Pecaut & Mamajek dwarf colour-temperature
/// table plus the solar anchor. It reproduces those points exactly and stays
/// within a few percent across O..M, which is ample for a *display* colour.
///
/// This is deliberately the simple first cut. The more accurate route is a
/// synthetic-passband "colour temperature" (invert the BP-RP a blackbody would
/// show through the real Gaia BP/RP filters); see `docs/star-colour.md`.
fn colour_temperature_from_bp_rp(bp_rp: f64) -> f64 {
    const A: f64 = 8235.0;
    const B: f64 = 0.997;
    const C: f64 = 1240.0;
    // Clamp the index to the range the fit was anchored over so a stray very
    // blue value can't drive the denominator toward zero (and T to infinity).
    let bp_rp = bp_rp.clamp(-0.5, 5.0);
    A / (bp_rp + B) + C
}

/// The XYZ radiance a single star contributes to the gather: its Pogson flux
/// scaled by its Y-normalised blackbody chromaticity, so `Y == flux` and X, Z
/// carry the hue. `redshift = 1.0` here (rest-frame); the black-hole redshift
/// is applied later at the tube via `temperature`.
fn star_emission_xyz(g_mag: f64, temperature: f64) -> CIETristimulus {
    let flux = pogson_flux(g_mag);
    let bb = get_cie_xyz_of_black_body_redshifted(temperature, 1.0);
    // Normalise out the absolute Planck amplitude (which for the disc encodes
    // brightness, but for a star must come from its magnitude, not its
    // temperature): keep only the chromaticity by dividing through by Y.
    let inv_y = if bb.y > 0.0 { 1.0 / bb.y } else { 0.0 };
    CIETristimulus::new(bb.x * inv_y * flux, flux, bb.z * inv_y * flux, 1.0)
}

pub struct StarCatalog {
    pub stars: Vec<Star>,
}

impl StarCatalog {
    /// Load a catalogue from a Parquet file on disk (e.g. `data/gaia_dr3.parquet`).
    pub fn load_parquet(path: impl AsRef<Path>) -> Result<Self, StarCatalogError> {
        Self::from_reader(File::open(path)?)
    }

    /// Load from any Parquet source (a file, or in-memory bytes for tests).
    pub fn from_reader<R: ChunkReader + 'static>(reader: R) -> Result<Self, StarCatalogError> {
        let batch_reader = ParquetRecordBatchReaderBuilder::try_new(reader)?
            .with_batch_size(8192)
            .build()?;

        let mut stars = Vec::new();
        for batch in batch_reader {
            let batch = batch?;
            // Columns are looked up by name (not position) so a change in query
            // column order can't silently mis-map them.
            let source_id = int64_column(&batch, "source_id")?;
            let ra = float64_column(&batch, "ra")?;
            let dec = float64_column(&batch, "dec")?;
            let g = float32_column(&batch, "phot_g_mean_mag")?;
            let bp = float32_column(&batch, "phot_bp_mean_mag")?;
            let rp = float32_column(&batch, "phot_rp_mean_mag")?;
            let bp_rp = float32_column(&batch, "bp_rp")?;

            stars.reserve(batch.num_rows());
            for i in 0..batch.num_rows() {
                let ra_deg = ra.value(i);
                let dec_deg = dec.value(i);
                let g_mag = g.value(i) as f64;
                let bp_rp_value = bp_rp.value(i) as f64;
                let temperature = colour_temperature_from_bp_rp(bp_rp_value);
                stars.push(Star {
                    source_id: source_id.value(i),
                    ra_deg,
                    dec_deg,
                    g_mag,
                    bp_mag: bp.value(i) as f64,
                    rp_mag: rp.value(i) as f64,
                    bp_rp: bp_rp_value,
                    direction: radec_to_unit_vector(ra_deg, dec_deg),
                    temperature,
                    emission_xyz: star_emission_xyz(g_mag, temperature),
                });
            }
        }

        Ok(Self { stars })
    }

    pub fn len(&self) -> usize {
        self.stars.len()
    }

    pub fn is_empty(&self) -> bool {
        self.stars.is_empty()
    }
}

/// ICRS (ra, dec) in degrees to a unit direction vector.
///
/// `dec` is a latitude (measured from the equator), while a polar angle `theta`
/// is a colatitude (measured from the pole): `theta = pi/2 - dec`. Substituting
/// that into the usual spherical-to-cartesian form gives the equatorial vector
/// directly, so this matches the renderer's escape directions.
pub fn radec_to_unit_vector(ra_deg: f64, dec_deg: f64) -> Vector3<f64> {
    let ra = ra_deg.to_radians();
    let dec = dec_deg.to_radians();
    Vector3::new(dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin())
}

fn column<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &'static str,
) -> Result<&'a arrow::array::ArrayRef, StarCatalogError> {
    batch
        .column_by_name(name)
        .ok_or(StarCatalogError::MissingColumn(name))
}

fn int64_column<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &'static str,
) -> Result<&'a Int64Array, StarCatalogError> {
    let col = column(batch, name)?;
    col.as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| StarCatalogError::UnexpectedColumnType(name, col.data_type().clone()))
}

fn float64_column<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &'static str,
) -> Result<&'a Float64Array, StarCatalogError> {
    let col = column(batch, name)?;
    col.as_any()
        .downcast_ref::<Float64Array>()
        .ok_or_else(|| StarCatalogError::UnexpectedColumnType(name, col.data_type().clone()))
}

fn float32_column<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &'static str,
) -> Result<&'a Float32Array, StarCatalogError> {
    let col = column(batch, name)?;
    col.as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| StarCatalogError::UnexpectedColumnType(name, col.data_type().clone()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn radec_maps_to_expected_axes() {
        // (ra=0, dec=0) points along +x; (ra=90, dec=0) along +y; dec=90 along +z.
        assert_relative_eq!(
            radec_to_unit_vector(0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            radec_to_unit_vector(90.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            radec_to_unit_vector(0.0, 90.0),
            Vector3::new(0.0, 0.0, 1.0),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            radec_to_unit_vector(0.0, -90.0),
            Vector3::new(0.0, 0.0, -1.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn radec_output_is_a_unit_vector() {
        for (ra, dec) in [(12.3, -45.6), (270.0, 12.0), (359.9, 89.9), (180.0, 0.0)] {
            assert_relative_eq!(radec_to_unit_vector(ra, dec).norm(), 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn relative_flux_is_pogson() {
        let star = Star {
            source_id: 1,
            ra_deg: 0.0,
            dec_deg: 0.0,
            g_mag: 5.0,
            bp_mag: 5.0,
            rp_mag: 5.0,
            bp_rp: 0.0,
            direction: Vector3::new(1.0, 0.0, 0.0),
            temperature: 5772.0,
            emission_xyz: CIETristimulus::new(0.0, 0.0, 0.0, 1.0),
        };
        // A 5-magnitude difference is exactly a factor of 100 in flux.
        assert_relative_eq!(star.relative_flux() * 100.0, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn colour_temperature_hits_reference_anchors() {
        // The reciprocal fit is pinned to these three (BP-RP, Teff) points.
        assert_relative_eq!(colour_temperature_from_bp_rp(0.0), 9500.0, epsilon = 5.0);
        assert_relative_eq!(colour_temperature_from_bp_rp(0.82), 5772.0, epsilon = 5.0);
        assert_relative_eq!(colour_temperature_from_bp_rp(3.0), 3300.0, epsilon = 5.0);
        // Monotone: bluer (smaller bp_rp) is hotter.
        assert!(colour_temperature_from_bp_rp(0.4) > colour_temperature_from_bp_rp(1.5));
    }

    #[test]
    fn emission_luminance_equals_flux_and_carries_hue() {
        // Y must equal the Pogson flux (brightness comes from the magnitude),
        // and a hot star must be bluer (Z/X larger) than a cool one.
        let hot = star_emission_xyz(0.0, 9000.0);
        let cool = star_emission_xyz(0.0, 3500.0);
        assert_relative_eq!(hot.y, pogson_flux(0.0), epsilon = 1e-12);
        assert_relative_eq!(cool.y, pogson_flux(0.0), epsilon = 1e-12);
        assert!((hot.z / hot.x) > (cool.z / cool.x));
    }

    /// Loads the real downloaded catalogue if present. Ignored by default
    /// because `data/` is git-ignored; run with `cargo test -- --ignored`.
    #[test]
    #[ignore]
    fn loads_downloaded_catalogue() {
        let catalog = StarCatalog::load_parquet("data/gaia_dr3.parquet").unwrap();
        assert!(!catalog.is_empty());
        for star in &catalog.stars {
            assert_relative_eq!(star.direction.norm(), 1.0, epsilon = 1e-9);
        }
    }
}
