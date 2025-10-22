use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Configuration for tone mapping operators.
///
/// Tone mapping compresses high dynamic range (HDR) values into low dynamic range (LDR)
/// for display, while preserving local contrast and detail.
///
/// This is essential for black hole rendering where brightness can vary by orders of
/// magnitude (e.g., bright accretion disk vs. dim background stars).
#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Debug)]
pub enum ToneMappingConfig {
    /// No tone mapping, just linear exposure multiplication.
    /// WARNING: Can cause clipping in high-contrast scenes.
    Linear { exposure: f64 },

    /// Reinhard tone mapping operator.
    ///
    /// Simple and fast, maps infinite input range to [0,1].
    /// Formula: L_out = L_in / (1 + L_in)
    ///
    /// Parameters:
    /// - exposure: Pre-exposure multiplier (typical: 0.5 - 2.0)
    /// - white_point: Luminance value that maps to white (optional, for local adaptation)
    ///
    /// Reference: Reinhard et al. "Photographic Tone Reproduction for Digital Images" (2002)
    Reinhard {
        exposure: f64,
        white_point: Option<f64>,
    },

    /// ACES filmic tone mapping (Academy Color Encoding System).
    ///
    /// Industry standard used in film production. Provides a film-like response curve
    /// with smooth highlight rolloff and good color preservation.
    ///
    /// This is a simplified approximation of the full ACES RRT (Reference Rendering Transform).
    ///
    /// Parameters:
    /// - exposure: Pre-exposure multiplier (typical: 1.0 - 3.0)
    ///
    /// Reference: Narkowicz "ACES Filmic Tone Mapping Curve" (2015)
    /// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    ACES { exposure: f64 },

    /// Uncharted 2 tone mapping (John Hable's filmic curve).
    ///
    /// Provides cinematic look with strong contrast and desaturated highlights.
    /// Used in the game Uncharted 2, good for dramatic scenes.
    ///
    /// Parameters:
    /// - exposure: Pre-exposure multiplier (typical: 1.0 - 4.0)
    ///
    /// Reference: Hable "Filmic Tonemapping Operators" (2010)
    /// http://filmicworlds.com/blog/filmic-tonemapping-operators/
    Uncharted2 { exposure: f64 },

    /// AgX tone mapping operator.
    ///
    /// Modern, perceptually-motivated tone mapper with excellent color preservation.
    /// Developed by Troy Sobotka as an open-source alternative to ACES.
    ///
    /// Provides more saturated colors than ACES while maintaining highlight rolloff.
    ///
    /// Parameters:
    /// - exposure: Pre-exposure multiplier (typical: 1.0 - 2.0)
    ///
    /// Reference: Sobotka "AgX" (2022)
    /// https://github.com/sobotka/AgX-S2O3
    AgX { exposure: f64 },
}

impl Default for ToneMappingConfig {
    fn default() -> Self {
        ToneMappingConfig::ACES { exposure: 1.0 }
    }
}

/// Applies tone mapping to linear RGB values.
///
/// Takes HDR linear RGB and returns tone-mapped linear RGB (before gamma correction).
/// The output should then be passed through sRGB gamma correction.
pub fn apply_tone_mapping(linear_rgb: Vector3<f64>, config: ToneMappingConfig) -> Vector3<f64> {
    match config {
        ToneMappingConfig::Linear { exposure } => linear_exposure(linear_rgb, exposure),
        ToneMappingConfig::Reinhard {
            exposure,
            white_point,
        } => reinhard(linear_rgb, exposure, white_point),
        ToneMappingConfig::ACES { exposure } => aces_filmic(linear_rgb, exposure),
        ToneMappingConfig::Uncharted2 { exposure } => uncharted2_filmic(linear_rgb, exposure),
        ToneMappingConfig::AgX { exposure } => agx(linear_rgb, exposure),
    }
}

/// Linear exposure: simply multiply by exposure value.
///
/// This is the simplest "tone mapping" but can cause severe clipping.
/// Only recommended for testing or very controlled scenes.
fn linear_exposure(rgb: Vector3<f64>, exposure: f64) -> Vector3<f64> {
    rgb * exposure
}

/// Reinhard tone mapping operator.
///
/// Global operator: L_out = L_in / (1 + L_in)
/// With white point: L_out = L_in * (1 + L_in/L_white²) / (1 + L_in)
///
/// Reference: Reinhard et al. SIGGRAPH 2002
fn reinhard(rgb: Vector3<f64>, exposure: f64, white_point: Option<f64>) -> Vector3<f64> {
    let color = rgb * exposure;

    match white_point {
        Some(white) => {
            // Extended Reinhard with white point
            // L_out = L * (1 + L/L_white²) / (1 + L)
            let white_sq = white * white;
            Vector3::new(
                reinhard_extended(color.x, white_sq),
                reinhard_extended(color.y, white_sq),
                reinhard_extended(color.z, white_sq),
            )
        }
        None => {
            // Simple Reinhard: L_out = L / (1 + L)
            Vector3::new(
                color.x / (1.0 + color.x),
                color.y / (1.0 + color.y),
                color.z / (1.0 + color.z),
            )
        }
    }
}

#[inline]
fn reinhard_extended(color: f64, white_sq: f64) -> f64 {
    color * (1.0 + color / white_sq) / (1.0 + color)
}

/// ACES filmic tone mapping curve (approximation).
///
/// Fits the full ACES RRT/ODT pipeline to a simple function.
/// Provides good contrast and color preservation.
///
/// Reference: Narkowicz 2015
/// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
fn aces_filmic(rgb: Vector3<f64>, exposure: f64) -> Vector3<f64> {
    let color = rgb * exposure;

    // ACES curve fitted parameters
    let a = 2.51;
    let b = 0.03;
    let c = 2.43;
    let d = 0.59;
    let e = 0.14;

    Vector3::new(
        aces_curve(color.x, a, b, c, d, e),
        aces_curve(color.y, a, b, c, d, e),
        aces_curve(color.z, a, b, c, d, e),
    )
}

#[inline]
fn aces_curve(x: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    ((x * (a * x + b)) / (x * (c * x + d) + e)).clamp(0.0, 1.0)
}

/// Uncharted 2 filmic tone mapping.
///
/// John Hable's filmic curve used in Uncharted 2.
/// Provides strong cinematic contrast with desaturated highlights.
///
/// Reference: Hable 2010
/// http://filmicworlds.com/blog/filmic-tonemapping-operators/
fn uncharted2_filmic(rgb: Vector3<f64>, exposure: f64) -> Vector3<f64> {
    let exposure_bias = 2.0; // Additional bias for proper white point

    let color = rgb * exposure;

    let mapped = Vector3::new(
        uncharted2_tonemap_impl(color.x * exposure_bias),
        uncharted2_tonemap_impl(color.y * exposure_bias),
        uncharted2_tonemap_impl(color.z * exposure_bias),
    );

    let white_scale = 1.0 / uncharted2_tonemap_impl(11.2);

    mapped * white_scale
}

#[inline]
fn uncharted2_tonemap_impl(x: f64) -> f64 {
    // Curve parameters (from Hable's presentation)
    let a = 0.15; // Shoulder strength
    let b = 0.50; // Linear strength
    let c = 0.10; // Linear angle
    let d = 0.20; // Toe strength
    let e = 0.02; // Toe numerator
    let f = 0.30; // Toe denominator

    ((x * (a * x + c * b) + d * e) / (x * (a * x + b) + d * f)) - e / f
}

/// AgX tone mapping operator.
///
/// This is a simplified version of Troy Sobotka's AgX.
/// The full AgX involves a more complex LUT-based transform.
/// This approximation provides similar characteristics.
///
/// Reference: https://github.com/sobotka/AgX-S2O3
fn agx(rgb: Vector3<f64>, exposure: f64) -> Vector3<f64> {
    let color = rgb * exposure;

    // AgX-like S-curve (simplified approximation)
    // Uses a smooth sigmoid-like function with good color preservation
    Vector3::new(
        agx_curve(color.x),
        agx_curve(color.y),
        agx_curve(color.z),
    )
}

#[inline]
fn agx_curve(x: f64) -> f64 {
    // Generalized logistic function (sigmoid variant)
    // Provides smooth rolloff similar to film response
    let x = x.max(0.0);

    // Parameters tuned for good highlight preservation
    let a = 1.0; // Maximum value
    let k = 2.0; // Steepness
    let b = 1.5; // Midpoint

    a / (1.0 + (x / b).powf(-k))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_linear_exposure_doubles_values() {
        let rgb = Vector3::new(0.5, 0.3, 0.2);
        let result = linear_exposure(rgb, 2.0);
        assert_abs_diff_eq!(result, Vector3::new(1.0, 0.6, 0.4));
    }

    #[test]
    fn test_reinhard_clamps_infinite_to_one() {
        // Reinhard should map very large values close to 1
        let rgb = Vector3::new(1000.0, 1000.0, 1000.0);
        let result = reinhard(rgb, 1.0, None);
        assert!(result.x > 0.99);
        assert!(result.x <= 1.0);
    }

    #[test]
    fn test_reinhard_preserves_blacks() {
        // Reinhard should preserve dark values
        let rgb = Vector3::new(0.01, 0.01, 0.01);
        let result = reinhard(rgb, 1.0, None);
        assert_abs_diff_eq!(result.x, 0.01 / 1.01, epsilon = 1e-6);
    }

    #[test]
    fn test_aces_in_range() {
        // ACES should always output values in [0, 1]
        let test_cases = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.5, 0.5, 0.5),
            Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(10.0, 10.0, 10.0),
            Vector3::new(100.0, 100.0, 100.0),
        ];

        for rgb in test_cases {
            let result = aces_filmic(rgb, 1.0);
            assert!(result.x >= 0.0 && result.x <= 1.0);
            assert!(result.y >= 0.0 && result.y <= 1.0);
            assert!(result.z >= 0.0 && result.z <= 1.0);
        }
    }

    #[test]
    fn test_uncharted2_in_range() {
        let rgb = Vector3::new(5.0, 3.0, 1.0);
        let result = uncharted2_filmic(rgb, 1.0);
        assert!(result.x >= 0.0 && result.x <= 1.0);
        assert!(result.y >= 0.0 && result.y <= 1.0);
        assert!(result.z >= 0.0 && result.z <= 1.0);
    }

    #[test]
    fn test_agx_monotonic() {
        // AgX should be monotonically increasing
        let values = vec![0.0, 0.5, 1.0, 2.0, 5.0];
        for i in 0..values.len() - 1 {
            let a = agx_curve(values[i]);
            let b = agx_curve(values[i + 1]);
            assert!(b >= a, "AgX should be monotonic: {} >= {}", b, a);
        }
    }

    #[test]
    fn test_tone_mapping_config_default() {
        let config = ToneMappingConfig::default();
        assert_eq!(config, ToneMappingConfig::ACES { exposure: 1.0 });
    }

    #[test]
    fn test_apply_tone_mapping_linear() {
        let rgb = Vector3::new(0.5, 0.3, 0.2);
        let config = ToneMappingConfig::Linear { exposure: 2.0 };
        let result = apply_tone_mapping(rgb, config);
        assert_abs_diff_eq!(result, Vector3::new(1.0, 0.6, 0.4));
    }

    #[test]
    fn test_apply_tone_mapping_aces() {
        let rgb = Vector3::new(1.0, 1.0, 1.0);
        let config = ToneMappingConfig::ACES { exposure: 1.0 };
        let result = apply_tone_mapping(rgb, config);
        // ACES should map 1.0 to somewhere around 0.75-0.85
        assert!(result.x > 0.7 && result.x < 0.9);
    }
}
