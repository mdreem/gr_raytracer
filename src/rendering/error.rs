//! Shared errors used throughout rendering.

use crate::rendering::camera::CameraError;
use crate::rendering::integration::integrator::IntegrationError;
use crate::rendering::texture::TextureError;
use image::ImageError;
use indicatif::style::TemplateError;
use std::io;

#[derive(Debug, thiserror::Error)]
pub enum RaytracerError {
    #[error("Integration error: {0}")]
    IntegrationError(#[from] IntegrationError),
    #[error("I/O error: {0}")]
    IoError(#[from] io::Error),
    #[error("Configuration file error: {0}")]
    ConfigurationFileError(io::Error),
    #[error("TOML error: {0}")]
    TomlError(#[from] toml::de::Error),
    #[error("Texture error: {0}")]
    TextureError(#[from] TextureError),
    #[error("Unphysical nonpositive redshift factor {0} in beaming")]
    UnphysicalRedshift(f64),
    #[error("Image error: {0}")]
    ImageError(#[from] ImageError),
    #[error("Image buffer creation failed")]
    ImageBufferCreation,
    #[error("Progress bar template error: {0}")]
    ProgressBarTemplateError(#[from] TemplateError),
    #[error("Camera error: {0}")]
    CameraError(#[from] CameraError),
    #[error("Invalid configuration: {0}")]
    InvalidConfiguration(String),
    #[error("No circular orbit possible")]
    NoCircularOrbitPossible,
    #[error("Radius is below RISCO")]
    BelowRISCO,
    #[error("Radius is not finite")]
    NonFiniteRadius,
    #[error("Number is below zero")]
    NumberBelowZero,
    #[error("Denominator is close to zero")]
    DenominatorCloseToZero,
    #[error("Invalid star tube: {0}")]
    InvalidStarTube(&'static str),
}
