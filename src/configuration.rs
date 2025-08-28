use serde::{Deserialize, Serialize};

#[derive(Deserialize, Serialize)]
pub struct RenderConfig {
    pub geometry_type: GeometryType,
    pub objects: Vec<ObjectsConfig>,
    pub celestial_texture: Option<String>,
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum GeometryType {
    Euclidean,
    EuclideanSpherical,
    Schwarzschild { radius: f64, horizon_epsilon: f64 },
}

#[derive(Deserialize, Serialize, Debug, PartialEq)]
pub enum TextureConfig {
    Bitmap {
        path: String,
    },
    Checker {
        width: f64,
        height: f64,
        color1: (u8, u8, u8),
        color2: (u8, u8, u8),
    },
    BlackBody {
        temperature: f64,
    },
}

#[derive(Deserialize, Serialize, Debug)]
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

    #[test]
    #[ignore]
    fn test_serialize() {
        let config = RenderConfig {
            celestial_texture: Some(String::from("resources/celestial_sphere.png")),
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
                    },
                },
                ObjectsConfig::Sphere {
                    radius: 0.5,
                    position: (4.4, 5.5, 6.6),
                    texture: TextureConfig::BlackBody {
                        temperature: 6500.0,
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
            celestial_texture = "resources/celestial_sphere.png"

            [geometry_type.Schwarzschild]
            radius = 2.0
            horizon_epsilon = 1e-4

            [[objects]]

            [objects.Sphere]
            radius = 1.0
            position = [1.1, 2.2, 3.3]

            [objects.Sphere.texture.Bitmap]
            path = "resources/sphere.png"

            [[objects]]

            [objects.Disc]
            inner_radius = 1.0
            outer_radius = 3.0

            [objects.Disc.texture.Checker]
            width = 0.5
            height = 0.5
            color1 = [255, 0, 0]
            color2 = [0, 0, 255]
        "#;

        let config: RenderConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(
            config.celestial_texture,
            Some(String::from("resources/celestial_sphere.png"))
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
                    path: String::from("resources/sphere.png")
                }
            );
        } else {
            panic!("Expected first object to be a Sphere");
        }

        if let ObjectsConfig::Disc {
            inner_radius,
            outer_radius,
            texture,
        } = &config.objects[1]
        {
            assert_eq!(*inner_radius, 1.0);
            assert_eq!(*outer_radius, 3.0);
        } else {
            panic!("Expected second object to be a Disc");
        }
    }
}
