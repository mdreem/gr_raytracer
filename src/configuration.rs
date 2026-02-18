use crate::rendering::color::CIETristimulusNormalization;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// Top-level scene configuration loaded from TOML.
#[derive(Deserialize, Serialize, JsonSchema, Clone)]
pub struct RenderConfig {
    /// Geometry model used for ray integration.
    pub geometry_type: GeometryType,
    /// Normalization applied to final rendered color values.
    pub color_normalization: CIETristimulusNormalization,
    /// Scene objects rendered in front of the celestial background.
    pub objects: Vec<ObjectsConfig>,
    /// Texture mapped onto the celestial sphere/background.
    pub celestial_texture: TextureConfig,
    /// Background emitter temperature in Kelvin.
    pub celestial_temperature: f64,
}

/// Available spacetime geometry models.
#[derive(Deserialize, Serialize, JsonSchema, Debug, PartialEq, Clone)]
pub enum GeometryType {
    /// Flat Euclidean geometry in Cartesian coordinates.
    Euclidean,
    /// Flat Euclidean geometry in spherical coordinates.
    EuclideanSpherical,
    /// Schwarzschild black hole geometry.
    Schwarzschild {
        /// Schwarzschild radius.
        radius: f64,
        /// Small offset used to avoid singular behavior at the horizon.
        horizon_epsilon: f64,
    },
    /// Kerr rotating black hole geometry.
    Kerr {
        /// Characteristic radius term.
        radius: f64,
        /// Kerr spin parameter.
        a: f64,
        /// Small offset used to avoid singular behavior at the horizon.
        horizon_epsilon: f64,
    },
}

/// Texture definitions used by objects and background.
#[derive(Deserialize, Serialize, JsonSchema, Debug, PartialEq, Clone)]
pub enum TextureConfig {
    /// Sample colors from a bitmap image.
    Bitmap {
        /// Exponent used for relativistic beaming intensity scaling.
        beaming_exponent: f64,
        /// Filesystem path to the image file.
        path: String,
        /// Normalization applied to sampled texture colors.
        color_normalization: CIETristimulusNormalization,
    },
    /// Procedural checkerboard texture.
    Checker {
        /// Exponent used for relativistic beaming intensity scaling.
        beaming_exponent: f64,
        /// Width of each checker cell in UV space.
        width: f64,
        /// Height of each checker cell in UV space.
        height: f64,
        /// RGB tuple for first checker color.
        color1: (u8, u8, u8),
        /// RGB tuple for second checker color.
        color2: (u8, u8, u8),
        /// Normalization applied to generated texture colors.
        color_normalization: CIETristimulusNormalization,
    },
    /// Color from black-body spectrum at emitter temperature.
    BlackBody {
        /// Exponent used for relativistic beaming intensity scaling.
        beaming_exponent: f64,
        /// Normalization applied to generated black-body colors.
        color_normalization: CIETristimulusNormalization,
    },
}

/// Scene object definitions.
#[derive(Deserialize, Serialize, JsonSchema, Debug, Clone)]
pub enum ObjectsConfig {
    /// Spherical emitter or textured sphere.
    Sphere {
        /// Sphere radius.
        radius: f64,
        /// Center position in Cartesian coordinates (x, y, z).
        position: (f64, f64, f64),
        /// Surface texture configuration.
        texture: TextureConfig,
        /// Emitter temperature in Kelvin.
        temperature: f64,
    },
    /// Thin disc emitter with inner and outer radius.
    Disc {
        /// Inner radius of the disc.
        inner_radius: f64,
        /// Outer radius of the disc.
        outer_radius: f64,
        /// Surface texture configuration.
        texture: TextureConfig,
        /// Emitter temperature in Kelvin.
        temperature: f64,
    },
    /// Participating-media disc with absorption/scattering.
    VolumetricDisc {
        /// Inner radius of the disc volume.
        inner_radius: f64,
        /// Outer radius of the disc volume.
        outer_radius: f64,
        /// Texture sampled within the volume.
        texture: TextureConfig,
        /// Base emitter temperature in Kelvin.
        temperature: f64,
        /// Optional rotation axis for the disc volume (x, y, z).
        axis: Option<(f64, f64, f64)>,
        /// Number of octaves used for procedural density noise.
        num_octaves: usize,
        /// Optional explicit seed for noise generation.
        perlin_seed: Option<u32>,
        /// Maximum integration steps through the participating media.
        max_steps: usize,
        /// Integration step size through the media volume.
        step_size: f64,
        /// Vertical thickness of the media disc.
        thickness: f64,
        /// Global multiplier applied to local density values.
        density_multiplier: f64,
        /// Reference temperature used for brightness scaling.
        brightness_reference_temperature: f64,
        /// Absorption coefficient.
        absorption: f64,
        /// Scattering coefficient.
        scattering: f64,
        /// Noise scaling factors (x, y, z).
        noise_scale: (f64, f64, f64),
        /// Scalar offset added to the sampled noise domain.
        noise_offset: f64,
    },
}

#[cfg(test)]
mod tests {
    use crate::configuration::{GeometryType, ObjectsConfig, RenderConfig, TextureConfig};
    use crate::rendering::color::CIETristimulusNormalization::NoNormalization;

    #[test]
    #[ignore]
    fn test_serialize() {
        let config = RenderConfig {
            celestial_texture: TextureConfig::Bitmap {
                beaming_exponent: 3.0,
                path: String::from("resources/celestial_sphere.png"),
                color_normalization: NoNormalization,
            },
            color_normalization: NoNormalization,
            celestial_temperature: 1500.0,
            geometry_type: GeometryType::Schwarzschild {
                radius: 2.0,
                horizon_epsilon: 1e-4,
            },
            objects: vec![
                ObjectsConfig::Sphere {
                    radius: 1.0,
                    position: (1.1, 2.2, 3.3),
                    texture: TextureConfig::Bitmap {
                        beaming_exponent: 3.0,
                        path: String::from("resources/sphere.png"),
                        color_normalization: NoNormalization,
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
                        color_normalization: NoNormalization,
                    },
                    temperature: 5500.0,
                },
                ObjectsConfig::Sphere {
                    radius: 0.5,
                    position: (4.4, 5.5, 6.6),
                    texture: TextureConfig::BlackBody {
                        beaming_exponent: 3.0,
                        color_normalization: NoNormalization,
                    },
                    temperature: 6500.0,
                },
            ],
        };

        let config_str = toml::to_string(&config).unwrap();
        println!("{}", config_str);
        assert!(false);
    }

    #[test]
    fn test_deserialize() {
        let toml_str = r#"
            color_normalization = "NoNormalization"
            celestial_temperature = 1500.0

            [celestial_texture.Bitmap]
            beaming_exponent = 3.0
            path = "resources/celestial_sphere.png"
            color_normalization = "NoNormalization"

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
            color_normalization = "NoNormalization"

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
            color_normalization = "NoNormalization"
        "#;

        let config: RenderConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(
            config.celestial_texture,
            TextureConfig::Bitmap {
                beaming_exponent: 3.0,
                path: String::from("resources/celestial_sphere.png"),
                color_normalization: NoNormalization,
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
                    color_normalization: NoNormalization,
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
