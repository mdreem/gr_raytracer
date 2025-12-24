use crate::rendering::raytracer::RaytracerError;
use log::{error, info};
use std::f64::consts::PI;

/// Trait for computing the temperature at a given radius.
pub trait TemperatureComputer: Sync {
    fn compute_temperature(&self, radius: f64) -> Result<f64, RaytracerError>;
}

/// A temperature computer that returns a constant temperature.
pub struct ConstantTemperatureComputer {
    temperature: f64,
}

impl ConstantTemperatureComputer {
    pub fn new(temperature: f64) -> Self {
        Self { temperature }
    }
}

impl TemperatureComputer for ConstantTemperatureComputer {
    fn compute_temperature(&self, _radius: f64) -> Result<f64, RaytracerError> {
        Ok(self.temperature)
    }
}

/// A temperature computer for a thin accretion disk around a Kerr black hole.
pub struct KerrTemperatureComputer {
    a: f64,
    radius: f64,
    /// Innermost stable circular orbit (ISCO) radius computed for the given
    /// `a` and `radius`.
    r_isco: f64,
    /// Mass accretion rate in internal units.
    m_dot: f64,
}

fn compute_r_isco(a: f64, radius: f64) -> f64 {
    let a_s = 2.0 * a / radius;

    let z1 = 1.0
        + (1.0 - a_s * a_s).powf(1.0 / 3.0)
            * ((1.0 + a_s).powf(1.0 / 3.0) + (1.0 - a_s).powf(1.0 / 3.0));
    let z2 = (3.0 * a_s * a_s + z1 * z1).sqrt();

    // prograde ISCO radius
    let r_isco = (3.0 + z2 - ((3.0 - z1) * (3.0 + z1 + 2.0 * z2)).sqrt()) * radius / 2.0;

    r_isco
}

const NUM_INTEGRATION_STEPS: i32 = 1000;
const NUM_STEPS_FIND_MAXIMUM: i32 = 10;

impl KerrTemperatureComputer {
    pub fn new(
        temperature: f64,
        outer_radius: f64,
        a: f64,
        radius: f64,
    ) -> Result<Self, RaytracerError> {
        let a_abs = a.abs(); // Ensure a co-rotating disc. TODO: generalize later.
        let r_isco = compute_r_isco(a_abs, radius);
        info!(
            "Computed r_isco: {} from a: {} and radius: {}",
            r_isco, a_abs, radius
        );
        let r0 = radius;

        let tmp_computer = Self {
            a: a_abs,
            radius,
            r_isco,
            m_dot: 1.0,
        };

        let mut max_f = 0.0;
        let mut max_r = 0.0;

        let dr = (outer_radius - r_isco) / NUM_STEPS_FIND_MAXIMUM as f64;
        for i in 0..NUM_STEPS_FIND_MAXIMUM {
            let r = r_isco + (i as f64 + 0.5) * dr;
            let f = tmp_computer.compute_f(r)?;
            if max_f < f {
                max_f = f;
                max_r = r;
            }
        }
        info!("Max f: {} at radius: {}", max_f, max_r);

        let integral = tmp_computer.compute_integral(max_r, r_isco)?;
        let pre_factor = tmp_computer.compute_prefactor(max_r)?;
        let coefficient = -1.0 / (PI * r0 * r0);

        // target: f = sigma_sb * T^4 = integral * pre_factor * coefficient * m_dot
        let sigma_sb = SIGMA_SB;
        let f = sigma_sb * temperature.powf(4.0);
        let m_dot = f / (coefficient * pre_factor * integral);
        info!(
            "Computed m_dot: {} for target temperature: {}",
            m_dot, temperature
        );

        Ok(Self {
            a: a_abs,
            radius,
            r_isco,
            m_dot,
        })
    }

    // TODO: Deduplicate against Kerr geometry methods.
    fn ut_contra(&self, r: f64) -> Result<f64, RaytracerError> {
        let a = self.a;
        let r_s = self.radius;
        let omega = self.angular_velocity(r);

        let g_tt = -(1.0 - r_s / r);
        let g_tphi = -a * r_s / r;
        let g_phiphi = r * r + a * a + a * a * r_s / r;

        let ut_pre = g_tt + 2.0 * omega * g_tphi + omega * omega * g_phiphi;

        if ut_pre >= 0.0 {
            error!(
                "No timelike circular orbit at r = {} (ut_pre = {})",
                r, ut_pre
            );
            return Err(RaytracerError::NoCircularOrbitPossible);
        }

        Ok((-ut_pre).sqrt().recip())
    }

    fn conserved_energy(&self, r: f64) -> Result<f64, RaytracerError> {
        let a = self.a;
        let r0 = self.radius;

        let g_tt = -1.0 + r0 / r;
        let g_tphi = -(r0 * a) / r;
        let omega = self.angular_velocity(r);

        let ut = self.ut_contra(r)?;

        Ok(-(g_tt + g_tphi * omega) * ut)
    }

    fn conserved_angular_momentum(&self, r: f64) -> Result<f64, RaytracerError> {
        let a = self.a;
        let r0 = self.radius;

        let g_tphi = -(r0 * a) / r;
        let g_phiphi = r * r + a * a + (r0 * a * a) / r;
        let omega = self.angular_velocity(r);
        let ut = self.ut_contra(r)?;

        Ok((g_tphi + g_phiphi * omega) * ut)
    }

    // TODO: Deduplicate against Kerr geometry methods.
    // https://arxiv.org/abs/1104.5499 equation (36)
    fn angular_velocity(&self, r: f64) -> f64 {
        let a = self.a;
        let r_s = self.radius;
        let m = 0.5 * r_s;
        let sqrt_m = m.sqrt();

        sqrt_m / (r.powf(1.5) + a * sqrt_m)
    }

    fn d_l_dr(&self, r: f64) -> Result<f64, RaytracerError> {
        let h = 1e-6 * r.max(1.0);

        if r - h < self.r_isco {
            // forward difference near ISCO to avoid going below ISCO
            Ok((self.conserved_angular_momentum(r + h)? - self.conserved_angular_momentum(r)?) / h)
        } else {
            Ok((self.conserved_angular_momentum(r + h)?
                - self.conserved_angular_momentum(r - h)?)
                / (2.0 * h))
        }
    }

    fn d_omega_dr(&self, r: f64) -> f64 {
        let h = 1e-10;
        let d_omega_dr = (self.angular_velocity(r + h) - self.angular_velocity(r - h)) / (2.0 * h);
        d_omega_dr
    }

    fn compute_f(&self, r: f64) -> Result<f64, RaytracerError> {
        let r_isco = self.r_isco;
        let r0 = self.radius;

        let integral = self.compute_integral(r, r_isco)?;
        let pre_factor = self.compute_prefactor(r)?;
        let coefficient = -self.m_dot / (PI * r0 * r0);

        Ok(coefficient * pre_factor * integral)
    }

    fn compute_prefactor(&self, r: f64) -> Result<f64, RaytracerError> {
        let e = self.conserved_energy(r)?;
        let l = self.conserved_angular_momentum(r)?;
        let omega = self.angular_velocity(r);
        let root_of_det_minus_g = r * r;

        let denominator = root_of_det_minus_g * (e - omega * l).powi(2);
        if denominator.abs() < 1e-20 {
            error!(
                "Denominator too small in pre_factor computation at r = {}: {}",
                r, denominator
            );
            return Err(RaytracerError::DenominatorCloseToZero);
        }
        let pre_factor = self.d_omega_dr(r) / denominator;
        Ok(pre_factor)
    }

    fn compute_integral(&self, r: f64, r_isco: f64) -> Result<f64, RaytracerError> {
        let dr = (r - r_isco) / NUM_INTEGRATION_STEPS as f64;
        let mut integral = 0.0;
        for i in 0..NUM_INTEGRATION_STEPS {
            let r_prime = r_isco + (i as f64 + 0.5) * dr;
            let e_prime = self.conserved_energy(r_prime)?;
            let l_prime = self.conserved_angular_momentum(r_prime)?;
            let omega_prime = self.angular_velocity(r_prime);

            integral += (e_prime - omega_prime * l_prime) * self.d_l_dr(r_prime)? * dr;
        }
        Ok(integral)
    }
}

const SIGMA_SB: f64 = 1.0; // Set to 1.0 as this will get calibrated later on anyway.

impl TemperatureComputer for KerrTemperatureComputer {
    fn compute_temperature(&self, radius: f64) -> Result<f64, RaytracerError> {
        if radius < self.r_isco {
            error!("Radius {} is below r_isco {}", radius, self.r_isco);
            return Err(RaytracerError::BelowRISCO);
        }

        let sigma_sb = SIGMA_SB;
        let f = self.compute_f(radius)?;
        if f < 0.0 {
            error!("Computed negative flux f: {} at radius: {}", f, radius);
            return Err(RaytracerError::NumberBelowZero);
        }
        let temp_fourth_power = f / sigma_sb;
        let temperature = temp_fourth_power.powf(0.25);
        Ok(temperature)
    }
}
