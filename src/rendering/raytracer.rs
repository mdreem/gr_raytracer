//! Render orchestration and focused sampling, star gathering, and output modules.

pub mod accuracy;
mod output;
mod progress;
mod sampling;
pub mod shadow_check;
mod star_gather;
pub(super) mod tube;

use crate::geometry::geometry::Geometry;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::integrator::StopReason;
use crate::rendering::light::color::ToneMappingMethod;
use crate::rendering::ray::IntegratedRay;
use crate::rendering::scene::Scene;
use log::{debug, info};

pub struct Raytracer<'a, G: Geometry> {
    pub scene: Scene<'a, G>,
    tone_mapping: ToneMappingMethod,
    exposure: f64,
}

impl<'a, G: Geometry> Raytracer<'a, G> {
    pub fn new(scene: Scene<'a, G>, tone_mapping: ToneMappingMethod, exposure: f64) -> Self {
        Self {
            scene,
            tone_mapping,
            exposure,
        }
    }

    #[allow(dead_code)] // For testing
    pub fn render_ray_at(&self, row: i64, col: i64) {
        let ray = self.scene.camera.get_ray_for(row, col);
        debug!("ray: {:?}", ray);
        let sample = self.scene.color_of_ray(&ray);
        debug!("sample: {:?}", sample);
    }

    pub fn integrate_ray_at_point(
        &self,
        row: i64,
        col: i64,
    ) -> Result<(IntegratedRay<G::Chart>, Option<StopReason>), RaytracerError> {
        let ray = self.scene.camera.get_ray_for(row, col);
        info!("ray for {}-{} is: {:?}", row, col, ray);
        self.scene.integrate_ray(&ray)
    }
}

#[cfg(test)]
mod tube_tests;

#[cfg(test)]
mod compositing_tests;

#[cfg(test)]
mod ring_test;
