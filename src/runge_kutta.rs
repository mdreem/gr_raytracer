use nalgebra::allocator::Allocator;
use nalgebra::OVector;
use nalgebra::{DefaultAllocator, Dim};

fn rk4<D: Dim>(
    x: &OVector<f64, D>,
    t: f64,
    h: f64,
    f: fn(f64, &OVector<f64, D>) -> OVector<f64, D>,
) -> OVector<f64, D>
where
    DefaultAllocator: Allocator<D>,
{
    let k1 = f(t, &x);
    let k2 = f(t + 0.5 * h, &(x + 0.5 * h * k1.clone()));
    let k3 = f(t + 0.5 * h, &(x + 0.5 * h * k2.clone()));
    let k4 = f(t + h, &(x + h * k3.clone()));

    x + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

#[cfg(test)]
mod tests {
    use crate::runge_kutta::rk4;
    use approx::assert_abs_diff_eq;
    use nalgebra::allocator::Allocator;
    use nalgebra::{Const, DefaultAllocator, OVector, Vector2};

    // d^2 y/dt^2 = a
    // -> dy_1/dt = y_2
    // -> dy_2/dt = a
    // initial conditions: y_1(0) = 1.0, y_2(0) = 2.0
    fn simple_equation(_: f64, y: &OVector<f64, Const<2>>) -> OVector<f64, Const<2>>
    where
        DefaultAllocator: Allocator<Const<2>>,
    {
        let a = 2.0;
        Vector2::new(y[1], a)
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

        let step_size = 0.0001;
        let mut t = 0.0;
        for i in 0..1000 {
            y = rk4(&y, t, step_size, simple_equation);
            t += step_size;
        }
        assert_abs_diff_eq!(y, solution_simple_equation(t), epsilon = 1e-6);

        for i in 0..1000 {
            y = rk4(&y, t, step_size, simple_equation);
            t += step_size;
        }
        assert_abs_diff_eq!(y, solution_simple_equation(t), epsilon = 1e-6);
    }
}
