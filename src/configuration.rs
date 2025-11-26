use crate::rendering::color::CIETristimulusNormalization;
use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize, Clone)]
pub struct RenderConfig {
    pub geometry_type: GeometryType,
    pub color_normalization: CIETristimulusNormalization,
    pub objects: Vec<ObjectsConfig>,
    pub celestial_texture: TextureConfig,
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Clone)]
pub enum GeometryType {
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
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Clone)]
pub enum TextureConfig {
    Bitmap {
        beaming_exponent: f64,
        path: String,
        color_normalization: CIETristimulusNormalization,
    },
    Checker {
        beaming_exponent: f64,
        width: f64,
        height: f64,
        color1: (u8, u8, u8),
        color2: (u8, u8, u8),
        color_normalization: CIETristimulusNormalization,
    },
    BlackBody {
        beaming_exponent: f64,
        temperature: f64,
        color_normalization: CIETristimulusNormalization,
    },
}

#[derive(Deserialize, Serialize, Debug, Clone)]
pub enum ObjectsConfig {
    Sphere {
        radius: f64,
        position: (f64, f64, f64),
        texture: TextureConfig,
    },
    Disc {
        inner_radius: f64,
        outer_radius: f64,
        texture: TextureConfig,
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
                path: String::from("resources/celestial_sphere.png"),
                color_normalization: NoNormalization,
            },
            color_normalization: NoNormalization,
            geometry_type: GeometryType::Schwarzschild {
                radius: 2.0,
                horizon_epsilon: 1e-4,
            },
            objects: vec![
                ObjectsConfig::Sphere {
                    radius: 1.0,
                    position: (1.1, 2.2, 3.3),
                    texture: TextureConfig::Bitmap {
                        path: String::from("resources/sphere.png"),
                        color_normalization: NoNormalization,
                    },
                },
                ObjectsConfig::Disc {
                    inner_radius: 1.0,
                    outer_radius: 3.0,
                    texture: TextureConfig::Checker {
                        width: 0.5,
                        height: 0.5,
                        color1: (255, 0, 0),
                        color2: (0, 0, 255),
                        color_normalization: NoNormalization,
                    },
                },
                ObjectsConfig::Sphere {
                    radius: 0.5,
                    position: (4.4, 5.5, 6.6),
                    texture: TextureConfig::BlackBody {
                        temperature: 6500.0,
                        color_normalization: NoNormalization,
                    },
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

            [celestial_texture.Bitmap]
            path = "resources/celestial_sphere.png"
            color_normalization = "NoNormalization"

            [geometry_type.Schwarzschild]
            radius = 2.0
            horizon_epsilon = 1e-4

            [[objects]]

            [objects.Sphere]
            radius = 1.0
            position = [1.1, 2.2, 3.3]

            [objects.Sphere.texture.Bitmap]
            path = "resources/sphere.png"
            color_normalization = "NoNormalization"

            [[objects]]

            [objects.Disc]
            inner_radius = 1.0
            outer_radius = 3.0

            [objects.Disc.texture.Checker]
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
                path: String::from("resources/celestial_sphere.png"),
                color_normalization: NoNormalization,
            }
        );

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
        } = &config.objects[0]
        {
            assert_eq!(*radius, 1.0);
            assert_eq!(position.0, 1.1);
            assert_eq!(position.1, 2.2);
            assert_eq!(position.2, 3.3);
            assert_eq!(
                texture,
                &TextureConfig::Bitmap {
                    path: String::from("resources/sphere.png"),
                    color_normalization: NoNormalization,
                }
            );
        } else {
            panic!("Expected first object to be a Sphere");
        }

        if let ObjectsConfig::Disc {
            inner_radius,
            outer_radius,
            texture: _,
        } = &config.objects[1]
        {
            assert_eq!(*inner_radius, 1.0);
            assert_eq!(*outer_radius, 3.0);
        } else {
            panic!("Expected second object to be a Disc");
        }
    }
}
