use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Default, Clone)]
pub struct RenderConfig {
    pub geometry_type: GeometryType,
    pub objects: Vec<ObjectsConfig>,
    pub celestial_texture: TextureConfig,
    pub celestial_temperature: f64,
    /// Which observer the camera is boosted to. Defaults to the static
    /// (Killing) observer so all geometries behave identically; use `Zamo`
    /// for close-in or ergosphere cameras (static observers do not exist
    /// there), or `Explicit` to pin a four-velocity in the geometry's native
    /// chart coordinates.
    #[serde(default)]
    pub camera_velocity: CameraVelocityConfig,
    /// Adaptive supersampling quality and edge-detection controls.
    #[serde(default)]
    pub adaptive_sampling: AdaptiveSamplingConfig,
    /// Optional Gaia DR3 star catalogue rendered as point sources on the
    /// celestial sphere. It is downloaded by scripts/gaia/download.py.
    /// When omitted, only `celestial_texture` is used.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub star_catalog: Option<StarCatalogConfig>,
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Clone)]
pub struct StarCatalogConfig {
    /// Path to the Parquet catalogue produced by `scripts/gaia/download.py`.
    pub path: String,
    /// Linear multiplier applied to each star's summed flux before it is
    /// written to the celestial sphere. Raise it to make the star field
    /// brighter relative to the tone-mapping exposure.
    #[serde(default = "default_star_flux_scale")]
    pub flux_scale: f64,
    /// Winding spread (in radians) across a traced tube's four corners above
    /// which the tube is subdivided before gathering stars. A tube whose
    /// corners have wound past this much relative angle straddles a caustic
    /// (photon ring), where the solid-angle magnification is unreliable, so it
    /// is split into four sub-tubes instead. Defaults to π.
    #[serde(default = "default_winding_spread_threshold")]
    pub winding_spread_threshold: f64,
    /// Maximum number of subdivision levels for a tube near the photon ring
    /// (0 disables subdivision, 1 permits one split). Each split traces five
    /// shared samples and creates four children; at most `4^depth` leaves.
    /// At the limit, escaped tubes use their finite area ratio and mixed
    /// tubes retain the corner-`a` fallback. Invalid areas or failed child
    /// rays are reported as render errors.
    #[serde(default = "default_max_subdivision_depth")]
    pub max_subdivision_depth: usize,
    /// Gather stars against a curved, arc-following tube boundary instead of
    /// the flat-quad corners, which fills gaps in the lensed star ring. Costs
    /// more near the critical curve, so it is opt-in for ring-heavy scenes. The
    /// `--curved-star-membership` CLI flag forces it on regardless of this.
    #[serde(default)]
    pub curved_star_membership: bool,
    /// Optional upper bound on a tube's lensing magnification `A/B`. Near a
    /// caustic (the photon/Einstein ring) the footprint collapses and the ratio
    /// blows up into resolution-dependent fireflies; a cap bounds it. Omit for
    /// no cap.
    #[serde(default)]
    pub magnification_cap: Option<f64>,
}

impl StarCatalogConfig {
    /// Reject parameter values that would silently corrupt the render: a
    /// non-finite or negative flux, a non-positive winding threshold, or a cap
    /// that is not strictly positive (a cap of 0 zeroes every tube's ratio and
    /// drops all stars).
    pub fn validate(&self) -> Result<(), String> {
        if !self.flux_scale.is_finite() || self.flux_scale < 0.0 {
            return Err(format!(
                "star_catalog.flux_scale must be finite and non-negative (got {})",
                self.flux_scale
            ));
        }
        if !self.winding_spread_threshold.is_finite() || self.winding_spread_threshold <= 0.0 {
            return Err(format!(
                "star_catalog.winding_spread_threshold must be finite and positive (got {})",
                self.winding_spread_threshold
            ));
        }
        if let Some(cap) = self.magnification_cap
            && (!cap.is_finite() || cap <= 0.0)
        {
            return Err(format!(
                "star_catalog.magnification_cap must be finite and positive (got {cap})"
            ));
        }
        Ok(())
    }
}

fn default_star_flux_scale() -> f64 {
    1.0
}

fn default_winding_spread_threshold() -> f64 {
    std::f64::consts::PI
}

fn default_max_subdivision_depth() -> usize {
    6
}

#[derive(Deserialize, Serialize, Clone, Debug, PartialEq)]
#[serde(default)]
pub struct AdaptiveSamplingConfig {
    pub enabled: bool,
    /// Number of strata along each pixel axis. The number of additional rays
    /// for a selected pixel is `samples_per_axis.pow(2)`.
    pub samples_per_axis: usize,
    pub luminance_contrast_threshold: f64,
    pub opacity_contrast_threshold: f64,
    /// Linear CIE Y floor below which contrast alone does not trigger sampling.
    /// `None` (the default) derives it per frame as a small fraction of the
    /// 99th-percentile luminance, so it tracks the disc brightness without a
    /// scene-specific constant. The derived value is section-relative, so a
    /// cropped section render is not byte-identical to the full frame there;
    /// set an explicit value for fully deterministic crops.
    pub minimum_luminance: Option<f64>,
    /// Accumulated object opacity required to classify a ray as an object hit.
    pub object_hit_opacity_threshold: f64,
    /// When true (default), pairs of `Escaped` (background / see-through) rays
    /// are never supersampled on contrast alone. This skips the high-frequency
    /// star field. Set false to also anti-alias a bright semi-transparent disc,
    /// at the cost of supersampling bright background stars too.
    pub exclude_background_contrast: bool,
}

impl Default for AdaptiveSamplingConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            samples_per_axis: 4,
            luminance_contrast_threshold: 0.15,
            opacity_contrast_threshold: 0.1,
            minimum_luminance: None,
            object_hit_opacity_threshold: 0.5,
            exclude_background_contrast: true,
        }
    }
}

impl AdaptiveSamplingConfig {
    pub fn validate(&self) -> Result<(), String> {
        if self.samples_per_axis == 0 {
            return Err("adaptive_sampling.samples_per_axis must be greater than zero".to_string());
        }
        for (name, value) in [
            (
                "luminance_contrast_threshold",
                self.luminance_contrast_threshold,
            ),
            (
                "opacity_contrast_threshold",
                self.opacity_contrast_threshold,
            ),
            (
                "object_hit_opacity_threshold",
                self.object_hit_opacity_threshold,
            ),
        ] {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(format!(
                    "adaptive_sampling.{name} must be finite and between 0 and 1 (got {value})"
                ));
            }
        }
        if let Some(value) = self.minimum_luminance
            && (!value.is_finite() || value < 0.0)
        {
            return Err(format!(
                "adaptive_sampling.minimum_luminance must be finite and non-negative (got {value})"
            ));
        }
        Ok(())
    }
}

#[derive(Deserialize, Serialize, Default, Debug, PartialEq, Clone)]
pub enum CameraVelocityConfig {
    /// Static (Killing) observer: at rest relative to infinity. Not defined
    /// inside the ergosphere of a spinning hole.
    #[default]
    StaticObserver,
    /// Zero angular momentum observer (locally non-rotating frame):
    /// co-rotates with frame dragging, exists everywhere outside the horizon.
    Zamo,
    /// Explicit four-velocity components in the geometry's native chart.
    /// Must be future-directed and normalized (u.u = +/-1 per signature).
    Explicit { components: [f64; 4] },
}

#[derive(Deserialize, Serialize, Default, Debug, PartialEq, Clone)]
pub enum GeometryType {
    #[default]
    Euclidean,
    EuclideanSpherical,
    Schwarzschild {
        radius: f64,
        horizon_epsilon: f64,
    },
    Kerr {
        radius: f64,
        a: f64,
        horizon_epsilon: f64,
    },
    KerrBL {
        radius: f64,
        a: f64,
        horizon_epsilon: f64,
    },
}

use crate::geometry::euclidean::EuclideanSpace;
use crate::geometry::euclidean_spherical::EuclideanSpaceSpherical;
use crate::geometry::geometry::RenderableGeometry;
use crate::geometry::kerr::Kerr;
use crate::geometry::kerr_bl::KerrBL;
use crate::geometry::schwarzschild::Schwarzschild;

impl GeometryType {
    pub fn get_renderable_geometry(&self) -> Box<dyn RenderableGeometry> {
        match self {
            GeometryType::Euclidean => Box::new(EuclideanSpace::new()),
            GeometryType::EuclideanSpherical => Box::new(EuclideanSpaceSpherical::new()),
            GeometryType::Schwarzschild {
                radius,
                horizon_epsilon,
            } => Box::new(Schwarzschild::new(*radius, *horizon_epsilon)),
            GeometryType::Kerr {
                radius,
                a,
                horizon_epsilon,
            } => Box::new(Kerr::new(*radius, *a, *horizon_epsilon)),
            GeometryType::KerrBL {
                radius,
                a,
                horizon_epsilon,
            } => Box::new(KerrBL::new(*radius, *a, *horizon_epsilon)),
        }
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Clone)]
pub enum TextureConfig {
    Bitmap {
        beaming_exponent: f64,
        path: String,
    },
    Checker {
        beaming_exponent: f64,
        width: f64,
        height: f64,
        color1: (u8, u8, u8),
        color2: (u8, u8, u8),
    },
    BlackBody {
        beaming_exponent: f64,
    },
}

impl Default for TextureConfig {
    fn default() -> Self {
        TextureConfig::BlackBody {
            beaming_exponent: 0.0,
        }
    }
}

fn default_disc_flux_scale() -> f64 {
    1.0
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub enum ObjectsConfig {
    Sphere {
        radius: f64,
        position: (f64, f64, f64),
        texture: TextureConfig,
        temperature: f64,
    },
    Disc {
        inner_radius: f64,
        outer_radius: f64,
        texture: TextureConfig,
        temperature: f64,
        #[serde(default = "default_disc_flux_scale")]
        flux_scale: f64,
    },
    VolumetricDisc {
        inner_radius: f64,
        outer_radius: f64,
        texture: TextureConfig,
        temperature: f64,
        axis: Option<(f64, f64, f64)>,
        num_octaves: usize,
        perlin_seed: Option<u32>,
        max_steps: usize,
        step_size: f64,
        thickness: f64,
        density_multiplier: f64,
        brightness_reference_temperature: f64,
        absorption: f64,
        scattering: f64,
        noise_scale: (f64, f64, f64),
        noise_offset: f64,
        #[serde(default = "default_disc_flux_scale")]
        flux_scale: f64,
    },
}

impl Default for ObjectsConfig {
    fn default() -> Self {
        ObjectsConfig::Sphere {
            radius: 1.0,
            position: (0.0, 0.0, 0.0),
            texture: TextureConfig::default(),
            temperature: 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::configuration::{
        AdaptiveSamplingConfig, GeometryType, ObjectsConfig, RenderConfig, TextureConfig,
    };

    #[test]
    fn adaptive_sampling_partial_config_uses_defaults() {
        let config: AdaptiveSamplingConfig = toml::from_str("samples_per_axis = 2").unwrap();
        assert_eq!(
            config,
            AdaptiveSamplingConfig {
                samples_per_axis: 2,
                ..Default::default()
            }
        );
        config.validate().unwrap();
    }

    #[test]
    fn adaptive_sampling_accepts_boundary_values() {
        AdaptiveSamplingConfig {
            luminance_contrast_threshold: 0.0,
            opacity_contrast_threshold: 1.0,
            minimum_luminance: Some(0.0),
            object_hit_opacity_threshold: 1.0,
            ..Default::default()
        }
        .validate()
        .unwrap();
    }

    #[test]
    fn adaptive_sampling_rejects_invalid_values() {
        let invalid_configs = [
            AdaptiveSamplingConfig {
                samples_per_axis: 0,
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                luminance_contrast_threshold: -0.1,
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                luminance_contrast_threshold: f64::NAN,
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                opacity_contrast_threshold: 1.1,
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                object_hit_opacity_threshold: f64::INFINITY,
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                minimum_luminance: Some(-0.1),
                ..Default::default()
            },
            AdaptiveSamplingConfig {
                minimum_luminance: Some(f64::INFINITY),
                ..Default::default()
            },
        ];

        for config in invalid_configs {
            assert!(config.validate().is_err(), "accepted {config:?}");
        }
    }

    #[test]
    fn test_serialize() {
        let config = RenderConfig {
            celestial_texture: TextureConfig::Bitmap {
                beaming_exponent: 3.0,
                path: String::from("resources/celestial_sphere.png"),
            },
            celestial_temperature: 1500.0,
            geometry_type: GeometryType::Schwarzschild {
                radius: 2.0,
                horizon_epsilon: 1e-4,
            },
            camera_velocity: Default::default(),
            adaptive_sampling: Default::default(),
            star_catalog: None,
            objects: vec![
                ObjectsConfig::Sphere {
                    radius: 1.0,
                    position: (1.1, 2.2, 3.3),
                    texture: TextureConfig::Bitmap {
                        beaming_exponent: 3.0,
                        path: String::from("resources/sphere.png"),
                    },
                    temperature: 4500.0,
                },
                ObjectsConfig::Disc {
                    inner_radius: 1.0,
                    outer_radius: 3.0,
                    texture: TextureConfig::Checker {
                        beaming_exponent: 3.0,
                        width: 0.5,
                        height: 0.5,
                        color1: (255, 0, 0),
                        color2: (0, 0, 255),
                    },
                    temperature: 5500.0,
                    flux_scale: 1.0,
                },
                ObjectsConfig::Sphere {
                    radius: 0.5,
                    position: (4.4, 5.5, 6.6),
                    texture: TextureConfig::BlackBody {
                        beaming_exponent: 3.0,
                    },
                    temperature: 6500.0,
                },
            ],
        };

        let config_str = toml::to_string(&config).expect("RenderConfig must be serialisable");

        // Round-trip via string: deserialising and re-serialising must produce
        // an identical TOML representation. RenderConfig and ObjectsConfig do
        // not derive PartialEq, so we compare the canonical TOML form rather
        // than the structs directly.
        let round_tripped: RenderConfig =
            toml::from_str(&config_str).expect("serialised RenderConfig must deserialise");
        let round_tripped_str =
            toml::to_string(&round_tripped).expect("round-tripped RenderConfig must reserialise");
        assert_eq!(config_str, round_tripped_str);
    }

    #[test]
    fn test_deserialize() {
        let toml_str = r#"
            celestial_temperature = 1500.0

            [celestial_texture.Bitmap]
            beaming_exponent = 3.0
            path = "resources/celestial_sphere.png"

            [geometry_type.Schwarzschild]
            radius = 2.0
            horizon_epsilon = 1e-4

            [[objects]]

            [objects.Sphere]
            radius = 1.0
            position = [1.1, 2.2, 3.3]
            temperature = 5500.0

            [objects.Sphere.texture.Bitmap]
            beaming_exponent = 3.0
            path = "resources/sphere.png"

            [[objects]]

            [objects.Disc]
            inner_radius = 1.0
            outer_radius = 3.0
            temperature = 6500.0

            [objects.Disc.texture.Checker]
            beaming_exponent = 3.0
            width = 0.5
            height = 0.5
            color1 = [255, 0, 0]
            color2 = [0, 0, 255]
        "#;

        let config: RenderConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(
            config.celestial_texture,
            TextureConfig::Bitmap {
                beaming_exponent: 3.0,
                path: String::from("resources/celestial_sphere.png"),
            }
        );
        assert_eq!(config.celestial_temperature, 1500.0);

        assert_eq!(
            config.geometry_type,
            GeometryType::Schwarzschild {
                radius: 2.0,
                horizon_epsilon: 1e-4,
            }
        );
        assert_eq!(config.objects.len(), 2);
        if let ObjectsConfig::Sphere {
            radius,
            position,
            texture,
            temperature,
        } = &config.objects[0]
        {
            assert_eq!(*radius, 1.0);
            assert_eq!(position.0, 1.1);
            assert_eq!(position.1, 2.2);
            assert_eq!(position.2, 3.3);
            assert_eq!(
                texture,
                &TextureConfig::Bitmap {
                    beaming_exponent: 3.0,
                    path: String::from("resources/sphere.png"),
                }
            );
            assert_eq!(*temperature, 5500.0);
        } else {
            panic!("Expected first object to be a Sphere");
        }

        if let ObjectsConfig::Disc {
            inner_radius,
            outer_radius,
            texture: _,
            temperature,
            ..
        } = &config.objects[1]
        {
            assert_eq!(*inner_radius, 1.0);
            assert_eq!(*outer_radius, 3.0);
            assert_eq!(*temperature, 6500.0);
        } else {
            panic!("Expected second object to be a Disc");
        }
    }
}
