use serde::{Deserialize, Serialize};

#[derive(Deserialize)]
pub struct RenderConfig {
    pub geometry_type: GeometryType,
    pub objects: Vec<ObjectsConfig>,
}

#[derive(Deserialize, Debug, PartialEq)]
pub enum GeometryType {
    Euclidean,
    EuclideanSpherical,
    Schwarzschild { radius: f64 },
}

#[derive(Deserialize, Debug)]
pub enum ObjectsConfig {
    Sphere {
        radius: f64,
    },
    Disc {
        inner_radius: f64,
        outer_radius: f64,
    },
}

#[cfg(test)]
mod tests {
    use crate::configuration::{GeometryType, ObjectsConfig, RenderConfig};

    #[test]
    fn test_deserialize() {
        let toml_str = r#"
            [geometry_type.Schwarzschild]
            radius = 2.0

            [[objects]]

            [objects.Sphere]
            radius = 1.0

            [[objects]]

            [objects.Disc]
            inner_radius = 1.0
            outer_radius = 3.0
        "#;

        let config: RenderConfig = toml::from_str(toml_str).unwrap();
        assert_eq!(
            config.geometry_type,
            GeometryType::Schwarzschild { radius: 2.0 }
        );
        assert_eq!(config.objects.len(), 2);
        if let ObjectsConfig::Sphere { radius } = &config.objects[0] {
            assert_eq!(*radius, 1.0);
        } else {
            panic!("Expected first object to be a Sphere");
        }

        if let ObjectsConfig::Disc {
            inner_radius,
            outer_radius,
        } = &config.objects[1]
        {
            assert_eq!(*inner_radius, 1.0);
            assert_eq!(*outer_radius, 3.0);
        } else {
            panic!("Expected second object to be a Disc");
        }
    }
}
