use crate::four_vector::CoordinateSystem::Spherical;
use crate::four_vector::{CoordinateSystem, FourVector};
use crate::geometry::{Geometry, HasCoordinateSystem, Tetrad};
use crate::runge_kutta::OdeFunction;
use crate::scene::EquationOfMotionState;
use nalgebra::{Const, Matrix4, OVector, Vector4};

#[derive(Clone, Debug)]
pub struct Schwarzschild {
    radius: f64,
}

impl Schwarzschild {
    pub fn new(radius: f64) -> Self {
        Schwarzschild { radius }
    }
}

impl OdeFunction<Const<8>> for Schwarzschild {
    // TODO: maybe just have geodesic being used in solver. This doesn't need to be that generic here.
    fn apply(&self, t: f64, y: &OVector<f64, Const<8>>) -> OVector<f64, Const<8>> {
        self.geodesic(t, y)
    }
}

impl HasCoordinateSystem for Schwarzschild {
    fn coordinate_system(&self) -> CoordinateSystem {
        Spherical
    }
}

// All coordinates here are spherical coordinates.
impl Geometry for Schwarzschild {
    fn geodesic(&self, _: f64, y: &EquationOfMotionState) -> EquationOfMotionState {
        let _t = y[0];
        let r = y[1];
        let theta = y[2];
        let _phi = y[3];

        let v_t = y[4];
        let v_r = y[5];
        let v_theta = y[6];
        let v_phi = y[7];

        let a = 1.0 - self.radius / r;
        let a_prime = self.radius / (r * r);
        let aprime_over_a = a_prime / a;

        // acceleration
        let a_t = -(aprime_over_a) * v_t * v_r;
        let a_r = -0.5 * a * a_prime * v_t * v_t
            + 0.5 * (aprime_over_a) * v_r * v_r
            + a * r * (v_theta * v_theta + v_phi * v_phi * theta.sin() * theta.sin());
        let a_theta = -(2.0 / r) * v_r * v_theta + theta.sin() * theta.cos() * v_phi * v_phi;
        let a_phi = -(2.0 / r) * v_phi * v_r - 2.0 * theta.cos() / theta.sin() * v_theta * v_phi;

        let y_new = EquationOfMotionState::from_column_slice(&[
            v_t, v_r, v_theta, v_phi, a_t, a_r, a_theta, a_phi,
        ]);
        y_new
    }

    // TODO: take into account rotations.
    fn get_tetrad_at(&self, position: &Vector4<f64>) -> Tetrad {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let rr0 = self.radius / r;
        let a = 1.0 - rr0;

        Tetrad::new(
            position.clone(),
            FourVector::new_spherical(1.0 / a, -rr0.sqrt(), 0.0, 0.0),
            FourVector::new_spherical(0.0, 0.0, 0.0, 1.0 / (r * theta.sin())), // Phi
            -FourVector::new_spherical(0.0, 0.0, 1.0 / r, 0.0),                // Theta
            -FourVector::new_spherical(-rr0.sqrt() / a, 1.0, 0.0, 0.0),        // R
        )
    }

    fn lorentz_transformation(
        &self,
        position: &Vector4<f64>,
        velocity: &FourVector,
    ) -> Matrix4<f64> {
        let mut matrix = Matrix4::zeros();
        let tetrad_t = self.get_tetrad_at(&position).t;

        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let a = 1.0 - self.radius / r;

        let metric_diag = Vector4::new(a, -1.0 / a, -r * r, -r * r * theta.sin() * theta.sin());

        let mut gamma = 0.0;
        for i in 0..4 {
            gamma += metric_diag[i] * velocity.vector[i] * tetrad_t.vector[i];
        }
        for mu in 0..4 {
            for nu in 0..4 {
                let mut res = 0.0;
                if mu == nu {
                    res = 1.0;
                }

                let a = 1.0 / (1.0 + gamma);
                let b = tetrad_t.vector[mu] + velocity.vector[mu];
                let c = metric_diag[nu] * (tetrad_t.vector[nu] + velocity.vector[nu]);
                res -= a * b * c;

                res += 2.0 * metric_diag[nu] * tetrad_t.vector[nu] * velocity.vector[mu];

                matrix[(mu, nu)] = res;
            }
        }
        matrix
    }

    fn mul(&self, position: &Vector4<f64>, v: &FourVector, w: &FourVector) -> f64 {
        let r = position[1];
        let theta = position[2];
        let _phi = position[3];

        let a = 1.0 - self.radius / r;

        a * v.vector[0] * w.vector[0]
            - v.vector[1] * w.vector[1] / a
            - r * r * v.vector[2] * w.vector[2]
            - r * r * theta.sin() * theta.sin() * v.vector[3] * w.vector[3]
    }
}

#[cfg(test)]
mod test_schwarzschild {
    use crate::four_vector::FourVector;
    use crate::runge_kutta::{rk4, OdeFunction};
    use crate::schwarzschild::Schwarzschild;
    use nalgebra::{Const, OVector, Vector2, Vector4};

    pub fn get_angular_momentum_from_phi(
        position: &Vector4<f64>,
        momentum: &FourVector,
        _geometry: &Schwarzschild,
    ) -> f64 {
        let r = position[1];
        let k_phi = momentum.vector[3];
        k_phi * r * r
    }

    pub fn get_energy_from_r(
        position: &Vector4<f64>,
        momentum: &FourVector,
        geometry: &Schwarzschild,
    ) -> f64 {
        let l = get_angular_momentum_from_phi(position, momentum, geometry);
        let r = position[1];
        let a = 1.0 - geometry.radius / r;
        let k_r = momentum.vector[1];
        let e_squared = k_r.powi(2) + a * l * l / (r * r);

        e_squared.sqrt()
    }

    pub fn get_energy_from_t(
        position: &Vector4<f64>,
        momentum: &FourVector,
        geometry: &Schwarzschild,
    ) -> f64 {
        let r = position[1];
        let a = 1.0 - geometry.radius / r;
        let k_t = momentum.vector[0];
        k_t * a
    }

    pub struct TestSchwarzschild {
        pub radius: f64,
    }

    #[derive(Debug)]
    pub struct Step {
        pub phi: f64,
        pub u: f64,
        du_dphi: f64,
    }

    impl OdeFunction<Const<2>> for TestSchwarzschild {
        fn apply(&self, _t: f64, y: &OVector<f64, Const<2>>) -> OVector<f64, Const<2>> {
            let u = y[0];
            let du_dphi = y[1];
            let d2u_dphi2 = -u + 3.0 * (self.radius / 2.0) * u.powi(2);
            Vector2::new(du_dphi, d2u_dphi2)
        }
    }

    impl TestSchwarzschild {
        pub fn compute_test_schwarzschild_trajectory(
            &self,
            r0: f64,
            l: f64,
            e: f64,
            max_steps: usize,
            step_size: f64,
        ) -> Vec<Step> {
            let u0 = 1.0 / r0;
            let b = l / e;
            // From energy equation before taking the second derivative.
            let du_dphi0 = (1.0 / b.powi(2) - u0.powi(2) * (1.0 - self.radius * u0)).sqrt();

            let mut y = Vector2::new(u0, du_dphi0);
            let mut t = 0.0;

            let mut result = vec![Step {
                phi: 0.0,
                u: u0,
                du_dphi: 0.0,
            }];

            for _ in 1..max_steps {
                y = rk4(&y, t, step_size, self);
                t += step_size;

                result.push(Step {
                    phi: t,
                    u: y[0],
                    du_dphi: y[1],
                });
            }

            result
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::camera::{Camera, Ray};
    use crate::debug::save_rays_to_file;
    use crate::four_vector::FourVector;
    use crate::geometry::Geometry;
    use crate::scene;
    use crate::scene::{CheckerMapper, Scene};
    use crate::schwarzschild::test_schwarzschild::{
        get_energy_from_r, get_energy_from_t, TestSchwarzschild,
    };
    use crate::schwarzschild::{test_schwarzschild, Schwarzschild};
    use crate::spherical_coordinates_helper::cartesian_to_spherical;
    use approx::assert_abs_diff_eq;
    use nalgebra::Vector4;
    use std::f64::consts::PI;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn test_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let geometry = Schwarzschild::new(2.0);

        let tetrad = geometry.get_tetrad_at(&position);

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.mul(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.t), 1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.x), -1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.y), -1.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.z, &tetrad.z), -1.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    #[test]
    fn test_lorentz_transformed_tetrad_orthonormal() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let r = position[1];
        let a = 1.0 - radius / position[1];

        let geometry = Schwarzschild::new(radius);

        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.

        let original_tetrad = geometry.get_tetrad_at(&position);
        let tetrad = crate::camera::lorentz_transform_tetrad(
            &geometry,
            &original_tetrad,
            &position,
            &velocity,
        );

        let k = tetrad.t + (-tetrad.z);
        let s = geometry.mul(&position, &k, &k);
        assert_abs_diff_eq!(s, 0.0);

        assert_abs_diff_eq!(
            tetrad.t.get_as_vector(),
            velocity.get_as_vector(),
            epsilon = 1e-6
        );

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.x), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.t, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.y), 0.0);
        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.x, &tetrad.z), 0.0);

        assert_abs_diff_eq!(geometry.mul(&position, &tetrad.y, &tetrad.z), 0.0);
    }

    #[ignore]
    #[test]
    fn save_camera_rays() {
        let rows = 30;
        let cols = 30;

        let position = cartesian_to_spherical(&Vector4::new(0.0, 0.0, 0.0, -10.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            rows,
            cols,
            Schwarzschild::new(2.0),
        );

        save_rays_to_file(rows, cols, &position, geometry, camera);
    }

    #[test]
    fn test_schwarzschild_ray() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 3.0, 4.0, 5.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            Schwarzschild::new(2.0),
        );

        let ray = camera.get_ray_for(1, 6);
        // This is a space-like vector and therefore should be -1.
        assert_abs_diff_eq!(
            geometry.mul(&position, &ray.direction, &ray.direction),
            -1.0
        );
        assert_abs_diff_eq!(geometry.mul(&position, &ray.momentum, &ray.momentum), 0.0);
    }

    fn create_camera(position: Vector4<f64>, radius: f64) -> Camera<Schwarzschild> {
        let r = position[1];
        let a = 1.0 - radius / r;
        let velocity = FourVector::new_spherical(1.0 / a, -(radius / r).sqrt(), 0.0, 0.0); // we have a freely falling observer here.
        let camera = Camera::new(
            position,
            velocity,
            PI / 2.0,
            11,
            11,
            Schwarzschild::new(radius),
        );
        camera
    }

    #[test]
    fn test_ray_null_condition_momentum() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 0.0, 0.0, 5.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let camera = create_camera(position, radius);

        for i in 1..11 {
            let ray = camera.get_ray_for(6, i);
            let m_s = geometry.mul(&position, &ray.momentum, &ray.momentum);
            assert_abs_diff_eq!(m_s, 0.0, epsilon = 1e-8);
        }
        for i in 1..11 {
            let ray = camera.get_ray_for(i, 6);
            let m_s = geometry.mul(&position, &ray.momentum, &ray.momentum);
            assert_abs_diff_eq!(m_s, 0.0, epsilon = 1e-8);
        }
    }

    #[test]
    fn test_ray_compare_conserved_quantities_in_equatorial_plane() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 0.0, 0.0, 5.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let camera = create_camera(position, radius);

        for i in 1..11 {
            let ray = camera.get_ray_for(6, i);
            assert_abs_diff_eq!(ray.direction.vector[2], 0.0); // ensure that the ray is in the equatorial plane.
            assert_abs_diff_eq!(ray.momentum.vector[2], 0.0); // ensure that the ray is in the equatorial plane.

            let m_s = geometry.mul(&position, &ray.momentum, &ray.momentum);
            assert_abs_diff_eq!(m_s, 0.0, epsilon = 1e-8);

            let e_r = get_energy_from_r(&position, &ray.momentum, &geometry);
            let e_t = get_energy_from_t(&position, &ray.momentum, &geometry);

            assert_abs_diff_eq!(e_r, e_t);
        }
    }

    #[test]
    fn test_trajectories_equal_with_rotated_momentum() {
        let position = cartesian_to_spherical(&Vector4::new(2.0, 0.0, 0.0, 5.0));
        let radius = 2.0;
        let geometry = Schwarzschild::new(radius);
        let camera = create_camera(position, radius);
        let scene: Scene<CheckerMapper, Schwarzschild> =
            scene::test_scene::create_scene(1.0, 2.0, 7.0, geometry.clone());

        let ray_a = camera.get_ray_for(6, 11);
        let ray_b = camera.get_ray_for(1, 6);

        // ensure rays are rotated by 90 degrees.
        assert_abs_diff_eq!(ray_a.momentum.vector[2], ray_b.momentum.vector[3]);
        assert_abs_diff_eq!(ray_a.momentum.vector[3], ray_b.momentum.vector[2]);

        println!("ray_a: {:?}", ray_a);
        println!("ray_b: {:?}", ray_b);

        let (trajectory_a, _) = scene.integrate(&ray_a);
        let (trajectory_b, _) = scene.integrate(&ray_b);
        assert_eq!(trajectory_a.len(), trajectory_b.len());

        for i in 0..trajectory_a.len() {
            let step_a = &trajectory_a[i];
            let step_b = &trajectory_b[i];

            assert_abs_diff_eq!(step_a.y[0], step_b.y[0], epsilon = 1e-8);
            assert_abs_diff_eq!(step_a.y[1], step_b.y[1], epsilon = 1e-8);
            assert_abs_diff_eq!(step_a.y[2], step_b.y[3], epsilon = 1e-8);
            assert_abs_diff_eq!(step_a.y[3], step_b.y[2], epsilon = 1e-8);
        }
    }

    #[test]
    fn compare_trajectories_escaping() {
        let radius = 1.0;

        let e = 1.0;
        let l = 5.0;

        let matching = compute_compared_trajectories(radius, e, l);

        assert_eq!(matching.len(), 204);
    }

    #[test]
    fn compare_trajectories_towards_black_hole() {
        let radius = 1.0;

        let e = 1.0;
        let l = 2.0;

        let matching = compute_compared_trajectories(radius, e, l);

        assert_eq!(matching.len(), 202);
    }

    fn compute_compared_trajectories(radius: f64, e: f64, l: f64) -> Vec<Point> {
        let position = Vector4::new(0.0, 5.0, PI / 2.0, 0.0);
        let r = position[1];
        let a = 1.0 - radius / r;

        let geometry = Schwarzschild::new(radius);
        let scene: Scene<CheckerMapper, Schwarzschild> =
            scene::test_scene::create_scene(1.0, 2.0, 7.0, geometry.clone());

        let momentum = FourVector::new_spherical(
            e / a,
            -(e * e - a * (l * l / (r * r))).sqrt(),
            0.0,
            l / (r * r),
        );

        let ray = Ray::new(0, 0, position, momentum, momentum);
        assert_abs_diff_eq!(geometry.mul(&position, &ray.momentum, &ray.momentum), 0.0);

        let ts = TestSchwarzschild { radius };

        let (result_geodesic_equation, _) = scene.integrate(&ray);
        let result_r_phi_equation = ts.compute_test_schwarzschild_trajectory(r, l, e, 204, 0.01);

        let matching = find_matching_points(
            &collect_points_step(&result_geodesic_equation),
            &collect_points_test_step(&result_r_phi_equation),
        );
        matching
    }

    fn save_trajectory(filename: &str, trajectory: &Vec<Point>) {
        let mut file = File::create(filename).expect("Unable to create file");
        file.write_all(b"r,phi\n").expect("Unable to write file");
        for step in trajectory {
            file.write_all(format!("{},{}\n", step.r, step.phi).as_bytes())
                .expect("Unable to write file");
        }
    }

    fn collect_points_step(steps: &Vec<scene::Step>) -> Vec<Point> {
        steps
            .iter()
            .map(|step| Point {
                r: step.y[1],
                phi: step.y[3],
            })
            .collect()
    }

    fn collect_points_test_step(steps: &Vec<test_schwarzschild::Step>) -> Vec<Point> {
        steps
            .iter()
            .map(|step| Point {
                r: step.u.recip(),
                phi: step.phi,
            })
            .collect()
    }

    fn find_matching_points(points_a: &Vec<Point>, points_b: &Vec<Point>) -> Vec<Point> {
        let mut result = Vec::new();

        let mut pos_b = 0;

        for point_a in points_a {
            for i in pos_b..points_b.len() {
                let point_b = &points_b[i];

                if (point_a.r - point_b.r).abs() < 0.1 && (point_a.phi - point_b.phi).abs() < 0.1 {
                    result.push(Point {
                        r: point_a.r,
                        phi: point_a.phi,
                    });
                    pos_b = i + 1;
                    break;
                }
            }
        }

        result
    }

    struct Point {
        pub r: f64,
        pub phi: f64,
    }
}
