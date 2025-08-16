use serde::Deserialize;

#[derive(Deserialize)]
pub struct RenderConfig {
    pub geometry_type: GeometryType,
    pub objects: Vec<ObjectsConfig>,
    pub celestial_texture: Option<String>,
}

#[derive(Deserialize, Debug, PartialEq)]
pub enum GeometryType {
    Euclidean,
    EuclideanSpherical,
    Schwarzschild { radius: f64, horizon_epsilon: f64 },
}

#[derive(Deserialize, Debug)]
pub enum ObjectsConfig {
    Sphere {
        radius: f64,
        position: (f64, f64, f64),
        texture: Option<String>,
    },
    Disc {
        inner_radius: f64,
        outer_radius: f64,
        texture: Option<String>,
    },
}

#[cfg(test)]
mod tests {
    use crate::configuration::{GeometryType, ObjectsConfig, RenderConfig};

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
            texture = "resources/sphere.png"

            [[objects]]

            [objects.Disc]
            inner_radius = 1.0
            outer_radius = 3.0
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
            assert_eq!(texture, &Some(String::from("resources/sphere.png")));
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
