//! Loader for the Gaia DR3 star-catalogue subset downloaded by
//! `scripts/gaia/download.py` (a Zstandard-compressed Parquet file).
//!
//! This is the renderer-side entry point of the point-source star pipeline
//! (see `docs/plan-12-star-catalog.md`). It reads the catalogue columns and
//! precomputes, per star, a unit direction on the celestial sphere. Working in
//! cartesian unit vectors keeps the later gather free of the (theta, phi)
//! coordinate singularities (pole and phi-seam): every membership / solid-angle
//! test is then a plain dot / cross product.

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
}

impl Star {
    /// Relative linear flux from the G magnitude (Pogson's law), normalized so
    /// magnitude 0 maps to flux 1. Absolute calibration (magnitude 0 -> chosen
    /// linear luminance at unit exposure) is a later concern of the gather.
    pub fn relative_flux(&self) -> f64 {
        10f64.powf(-0.4 * self.g_mag)
    }
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
                stars.push(Star {
                    source_id: source_id.value(i),
                    ra_deg,
                    dec_deg,
                    g_mag: g.value(i) as f64,
                    bp_mag: bp.value(i) as f64,
                    rp_mag: rp.value(i) as f64,
                    bp_rp: bp_rp.value(i) as f64,
                    direction: radec_to_unit_vector(ra_deg, dec_deg),
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
        assert_relative_eq!(radec_to_unit_vector(0.0, 0.0), Vector3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(radec_to_unit_vector(90.0, 0.0), Vector3::new(0.0, 1.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(radec_to_unit_vector(0.0, 90.0), Vector3::new(0.0, 0.0, 1.0), epsilon = 1e-12);
        assert_relative_eq!(radec_to_unit_vector(0.0, -90.0), Vector3::new(0.0, 0.0, -1.0), epsilon = 1e-12);
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
        };
        // A 5-magnitude difference is exactly a factor of 100 in flux.
        assert_relative_eq!(star.relative_flux() * 100.0, 1.0, epsilon = 1e-12);
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
