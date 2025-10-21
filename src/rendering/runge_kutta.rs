use crate::rendering::integrator::IntegrationError;
use crate::rendering::raytracer::RaytracerError;
use nalgebra::OVector;
use nalgebra::allocator::Allocator;
use nalgebra::{DefaultAllocator, Dim};

pub trait OdeFunction<D: Dim>
where
    DefaultAllocator: Allocator<D>,
{
    fn apply(&self, t: f64, y: &OVector<f64, D>) -> OVector<f64, D>;
}

// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method

const A1: f64 = 0.0;
const A2: f64 = 2.0 / 9.0;
const A3: f64 = 1.0 / 3.0;
const A4: f64 = 3.0 / 4.0;
const A5: f64 = 1.0;
const A6: f64 = 5.0 / 6.0;

const B_2_1: f64 = 2.0 / 9.0;
const B_3_1: f64 = 1.0 / 12.0;
const B_4_1: f64 = 69.0 / 128.0;
const B_5_1: f64 = -17.0 / 12.0;
const B_6_1: f64 = 65.0 / 432.0;

const B_3_2: f64 = 1.0 / 4.0;
const B_4_2: f64 = -243.0 / 128.0;
const B_5_2: f64 = 27.0 / 4.0;
const B_6_2: f64 = -5.0 / 16.0;

const B_4_3: f64 = 135.0 / 64.0;
const B_5_3: f64 = -27.0 / 5.0;
const B_6_3: f64 = 13.0 / 16.0;

const B_5_4: f64 = 16.0 / 15.0;
const B_6_4: f64 = 4.0 / 27.0;

const B_6_5: f64 = 5.0 / 144.0;

const C_1: f64 = 1.0 / 9.0;
const C_2: f64 = 0.0;
const C_3: f64 = 9.0 / 20.0;
const C_4: f64 = 16.0 / 45.0;
const C_5: f64 = 1.0 / 12.0;

const CH_1: f64 = 47.0 / 450.0;
const CH_2: f64 = 0.0;
const CH_3: f64 = 12.0 / 25.0;
const CH_4: f64 = 32.0 / 225.0;
const CH_5: f64 = 1.0 / 30.0;
const CH_6: f64 = 6.0 / 25.0;

const CT_1: f64 = 1.0 / 150.0;
const CT_2: f64 = 0.0;
const CT_3: f64 = -3.0 / 100.0;
const CT_4: f64 = 16.0 / 75.0;
const CT_5: f64 = 1.0 / 20.0;
const CT_6: f64 = -6.0 / 25.0;

const BETA: f64 = 0.9;
const CONVERGENCY_ORDER: f64 = 5.0;
const ERROR_RATIO_SMALL_ERROR: f64 = 1e-5;
const MAX_RETRY_STEP: i32 = 100;

fn rkf45_step<D: Dim>(
    y: &OVector<f64, D>,
    t: f64,
    h: f64,
    f: &dyn OdeFunction<D>,
) -> (OVector<f64, D>, f64)
where
    DefaultAllocator: Allocator<D>,
{
    let k1 = h * f.apply(t + A1 * h, &y);
    let k2 = h * f.apply(t + A2 * h, &(y + B_2_1 * k1.clone()));
    let k3 = h * f.apply(t + A3 * h, &(y + B_3_1 * k1.clone() + B_3_2 * k2.clone()));
    let k4 = h * f.apply(
        t + A4 * h,
        &(y + B_4_1 * k1.clone() + B_4_2 * k2.clone() + B_4_3 * k3.clone()),
    );
    let k5 = h * f.apply(
        t + A5 * h,
        &(y + B_5_1 * k1.clone() + B_5_2 * k2.clone() + B_5_3 * k3.clone() + B_5_4 * k4.clone()),
    );
    let k6 = h * f.apply(
        t + A6 * h,
        &(y + B_6_1 * k1.clone()
            + B_6_2 * k2.clone()
            + B_6_3 * k3.clone()
            + B_6_4 * k4.clone()
            + B_6_5 * k5.clone()),
    );

    let y_new = y
        + CH_1 * k1.clone()
        + CH_2 * k2.clone()
        + CH_3 * k3.clone()
        + CH_4 * k4.clone()
        + CH_5 * k5.clone()
        + CH_6 * k6.clone();
    let truncation_error =
        (CT_1 * k1 + CT_2 * k2 + CT_3 * k3 + CT_4 * k4 + CT_5 * k5 + CT_6 * k6).norm();
    (y_new, truncation_error)
}

pub fn rkf45<D: Dim>(
    y: &OVector<f64, D>,
    t: f64,
    h: f64,
    epsilon: f64,
    f: &dyn OdeFunction<D>,
) -> Result<(OVector<f64, D>, f64), RaytracerError>
where
    DefaultAllocator: Allocator<D>,
{
    let mut h_cur = h;
    for _i in 0..MAX_RETRY_STEP {
        let (y_new, truncation_error) = rkf45_step(y, t, h_cur, f);

        let h_new = BETA * h_cur * (epsilon / truncation_error).powf(1.0 / CONVERGENCY_ORDER);

        if truncation_error > epsilon {
            h_cur = h_new / 2.0; // Halve the step size and retry.
        } else {
            // step is accepted
            return if truncation_error / epsilon < ERROR_RATIO_SMALL_ERROR {
                Ok((y_new, 2.0 * h_cur))
            } else {
                Ok((y_new, h_cur))
            };
        }
    }
    Err(RaytracerError::IntegrationError(
        IntegrationError::MaxStepsReached,
    ))
}

#[cfg(test)]
mod tests {
    use crate::rendering::runge_kutta::{OdeFunction, rkf45};
    use approx::assert_abs_diff_eq;
    use nalgebra::{Const, OVector, Vector2};

    struct SimpleEquation {}

    // d^2 y/dt^2 = a
    // -> dy_1/dt = y_2
    // -> dy_2/dt = a
    // initial conditions: y_1(0) = 1.0, y_2(0) = 2.0
    impl OdeFunction<Const<2>> for SimpleEquation {
        fn apply(&self, _: f64, y: &OVector<f64, Const<2>>) -> OVector<f64, Const<2>> {
            let a = 2.0;
            Vector2::new(y[1], a)
        }
    }

    fn solution_simple_equation(t: f64) -> OVector<f64, Const<2>> {
        let a = 2.0;
        let x1_0 = 1.0;
        let x2_0 = 2.0;

        let x1 = 0.5 * a * t * t + x2_0 * t + x1_0;
        let x2 = a * t + x2_0;

        Vector2::new(x1, x2)
    }

    #[test]
    fn rk4_works() {
        let mut y = Vector2::new(1.0, 2.0);
        let simple_equation = SimpleEquation {};

        let step_size = 0.0000001;
        let mut t = 0.0;
        let mut h = step_size;
        while t <= 25.0 {
            (y, h) = rkf45(&y, t, h, 1e-10, &simple_equation).unwrap();
            t += h;
        }
        assert_abs_diff_eq!(y, solution_simple_equation(t - h), epsilon = 1e-5);

        while t <= 50.0 {
            (y, h) = rkf45(&y, t, h, 1e-10, &simple_equation).unwrap();
            t += h;
        }
        assert_abs_diff_eq!(y, solution_simple_equation(t - h), epsilon = 1e-5);
    }
}
