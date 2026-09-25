use crate::rendering::color::CIETristimulus;
use nalgebra::Vector3;

/// Integrated linear XYZ light and the fraction of background light transmitted.
///
/// XYZ is already weighted by opacity/coverage (premultiplied), unlike the
/// straight colors returned by textures. It must never be multiplied by this
/// layer's opacity again. HDR light is unbounded; transmittance is in [0, 1].
#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Radiance {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub transmittance: f64,
}

impl Radiance {
    pub const TRANSPARENT: Self = Self::new(0.0, 0.0, 0.0, 1.0);
    pub const BLACK: Self = Self::new(0.0, 0.0, 0.0, 0.0);

    /// Construct already-integrated light; no opacity multiplication occurs.
    pub const fn new(x: f64, y: f64, z: f64, transmittance: f64) -> Self {
        Self {
            x,
            y,
            z,
            transmittance,
        }
    }

    /// Convert a straight-alpha surface/sky texture at the material boundary.
    pub fn from_straight(color: CIETristimulus) -> Self {
        let opacity = color.alpha.clamp(0.0, 1.0);
        Self::new(
            color.x * opacity,
            color.y * opacity,
            color.z * opacity,
            1.0 - opacity,
        )
    }

    pub fn opacity(self) -> f64 {
        1.0 - self.transmittance
    }

    pub fn as_vector(self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Compose this foreground over a farther layer: L = Lf + Tf Lb, T = Tf Tb.
    pub fn over(self, background: Self) -> Self {
        Self::new(
            self.x + self.transmittance * background.x,
            self.y + self.transmittance * background.y,
            self.z + self.transmittance * background.z,
            self.transmittance * background.transmittance,
        )
    }

    /// Flatten onto black for RGB/HDR export. No unpremultiplication: the XYZ
    /// values are exactly the light reaching the pixel, even for transparent rays.
    pub fn to_xyz_over_black(self) -> CIETristimulus {
        CIETristimulus::new(self.x, self.y, self.z, 1.0)
    }

    /// Return a filtered surface texture to its straight-alpha interface.
    /// This is not used for traced radiance or final image output.
    pub fn to_straight_surface(self) -> CIETristimulus {
        let opacity = self.opacity();
        if opacity <= 0.0 {
            CIETristimulus::new(0.0, 0.0, 0.0, 0.0)
        } else {
            CIETristimulus::new(
                self.x / opacity,
                self.y / opacity,
                self.z / opacity,
                opacity,
            )
        }
    }
}

/// Area/sample averaging is distinct from front-to-back layer composition.
/// Average integrated light and transmittance, never straight RGB/XYZ.
#[derive(Default)]
pub struct RadianceMean {
    light: Vector3<f64>,
    transmittance: f64,
    weight: f64,
}

impl RadianceMean {
    pub fn add(&mut self, sample: Radiance, weight: f64) {
        self.light += sample.as_vector() * weight;
        self.transmittance += sample.transmittance * weight;
        self.weight += weight;
    }

    pub fn mean(self) -> Option<Radiance> {
        (self.weight > 0.0).then(|| {
            let light = self.light / self.weight;
            Radiance::new(light.x, light.y, light.z, self.transmittance / self.weight)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn slab(source: f64, tau: f64) -> Radiance {
        let light = source * -(-tau).exp_m1();
        Radiance::new(light, light, light, (-tau).exp())
    }

    #[test]
    fn straight_surface_is_premultiplied_once() {
        let foreground = Radiance::from_straight(CIETristimulus::new(0.6, 0.4, 0.2, 0.5));
        let background = Radiance::new(0.2, 0.4, 0.6, 0.0);
        let result = foreground.over(background);
        assert_abs_diff_eq!(result.as_vector(), Vector3::repeat(0.4));
        assert_eq!(result.transmittance, 0.0);
        assert_eq!(foreground.over(Radiance::TRANSPARENT), foreground);
    }

    #[test]
    fn transparent_surface_colors_do_not_contribute_light() {
        let hidden = Radiance::from_straight(CIETristimulus::new(100.0, 20.0, 30.0, 0.0));
        assert_eq!(hidden, Radiance::TRANSPARENT);
        assert_eq!(hidden.over(hidden), Radiance::TRANSPARENT);
        let background = Radiance::new(2.0, 3.0, 4.0, 0.0);
        assert_eq!(hidden.over(background), background);
    }

    #[test]
    fn slab_matches_analytic_transfer_over_black_and_bright_backgrounds() {
        for tau in [0.0, 1e-10, 0.1, 1.0, 20.0] {
            let foreground = slab(2.0, tau);
            for background in [0.0, 3.0, 1000.0] {
                let result =
                    foreground.over(Radiance::new(background, background, background, 0.0));
                let expected = 2.0 * -(-tau).exp_m1() + background * (-tau).exp();
                assert_abs_diff_eq!(
                    result.as_vector(),
                    Vector3::repeat(expected),
                    epsilon = 1e-12
                );
            }
        }
    }

    #[test]
    fn slab_splitting_and_composition_order_are_consistent() {
        let whole = slab(2.0, 0.8);
        let split = slab(2.0, 0.3).over(slab(2.0, 0.5));
        assert_abs_diff_eq!(whole.as_vector(), split.as_vector(), epsilon = 1e-14);
        assert_abs_diff_eq!(whole.transmittance, split.transmittance, epsilon = 1e-14);
        let near = slab(2.0, 0.3);
        let far = slab(5.0, 0.5);
        let background = Radiance::new(10.0, 10.0, 10.0, 0.0);
        assert_abs_diff_eq!(
            near.over(far).over(background).as_vector(),
            near.over(far.over(background)).as_vector(),
            epsilon = 1e-14
        );
        assert_abs_diff_eq!(
            near.over(far).y,
            2.0 * (1.0 - (-0.3_f64).exp()) + (-0.3_f64).exp() * 5.0 * (1.0 - (-0.5_f64).exp()),
            epsilon = 1e-14
        );
    }

    #[test]
    fn sample_average_preserves_coverage_and_background_light() {
        let covered = Radiance::from_straight(CIETristimulus::new(2.0, 1.0, 0.5, 1.0));
        let sky = Radiance::TRANSPARENT;
        let background = Radiance::new(4.0, 6.0, 8.0, 0.0);
        let mut average = RadianceMean::default();
        average.add(covered, 1.0);
        average.add(sky, 1.0);
        let average = average.mean().unwrap();
        assert_eq!(average.transmittance, 0.5);
        assert_abs_diff_eq!(average.as_vector(), Vector3::new(1.0, 0.5, 0.25));
        assert_abs_diff_eq!(
            average.over(background).as_vector(),
            0.5 * (covered.over(background).as_vector() + sky.over(background).as_vector())
        );
        assert_eq!(RadianceMean::default().mean(), None);
    }

    #[test]
    fn export_keeps_integrated_light_without_opacity_scaling_or_division() {
        let volume = Radiance::new(0.1, 0.2, 0.3, 0.9);
        let output = volume.to_xyz_over_black();
        assert_eq!(output.as_vector(), volume.as_vector());
        assert_eq!(output.alpha, 1.0);
    }
}
