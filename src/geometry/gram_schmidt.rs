use crate::geometry::four_vector::FourVector;
use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use std::ops::Mul;

fn proj<G: Geometry>(geometry: &G, position: &Point, u: &FourVector, v: &FourVector) -> FourVector {
    let p1 = geometry.inner_product(position, v, u);
    let p2 = geometry.inner_product(position, u, u);

    (p1 / p2) * u.clone()
}

pub fn gram_schmidt<G: Geometry>(
    geometry: &G,
    position: &Point,
    vectors: &[FourVector],
) -> Vec<FourVector> {
    let mut orthonormal_vectors: Vec<FourVector> = Vec::new();

    for v in vectors {
        let mut w = v.clone();

        for u in &orthonormal_vectors {
            let projection = proj(geometry, position, u, &w);
            w = w + (-projection);
        }

        let norm = geometry.inner_product(position, &w, &w).abs().sqrt();
        println!("norm = {:?}", norm);
        orthonormal_vectors.push((1.0 / norm) * w);
    }

    orthonormal_vectors
}

#[cfg(test)]
mod tests {
    use crate::geometry::euclidean::EuclideanSpace;
    use crate::geometry::four_vector::FourVector;
    use crate::geometry::geometry::InnerProduct;
    use crate::geometry::gram_schmidt::gram_schmidt;
    use crate::geometry::point::Point;
    use crate::geometry::schwarzschild::Schwarzschild;
    use crate::geometry::spherical_coordinates_helper::cartesian_to_spherical;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_gram_schmidt_properties_euclidean() {
        let position = Point::new_cartesian(0.0, 3.0, 2.0, 1.0);
        let point_in_spherical = position.get_as_spherical();
        let _r = point_in_spherical.x;
        let theta = point_in_spherical.y;
        let phi = point_in_spherical.z;

        let e_t = FourVector::new_cartesian(1.0, 0.0, 0.0, 0.0);
        let e_r = FourVector::new_cartesian(
            0.0,
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            theta.cos(),
        );
        let e_theta = FourVector::new_cartesian(0.0, 0.0, -1.0, 0.0);
        let e_phi = FourVector::new_cartesian(0.0, 0.0, 0.0, -1.0);

        let geometry = EuclideanSpace::new();
        let vectors = vec![e_t, e_r, e_theta, e_phi];
        let orthonormal_vectors = gram_schmidt(&geometry, &position, &vectors);

        println!("{:?}", orthonormal_vectors);

        let norm =
            geometry.inner_product(&position, &orthonormal_vectors[0], &orthonormal_vectors[0]);
        assert_abs_diff_eq!(norm, 1.0);
        assert_eq!(orthonormal_vectors.len(), 4);
        for i in 1..4 {
            let norm =
                geometry.inner_product(&position, &orthonormal_vectors[i], &orthonormal_vectors[i]);
            assert_abs_diff_eq!(norm, -1.0);
        }

        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    assert_abs_diff_eq!(
                        geometry.inner_product(
                            &position,
                            &orthonormal_vectors[i],
                            &orthonormal_vectors[j]
                        ),
                        0.0
                    );
                }
            }
        }
    }

    #[test]
    fn test_gram_schmidt_properties_schwarzschild() {
        let r_s = 2.0;

        let position = Point::new_cartesian(0.0, 3.0, 2.0, 1.0);
        let point_in_spherical = cartesian_to_spherical(&position);
        println!("{:?}", point_in_spherical);
        let r = point_in_spherical.vector[1];
        let theta = point_in_spherical.vector[2];
        let phi = point_in_spherical.vector[3];

        let a = 1.0 - r_s / r;

        let e_t = FourVector::new_spherical(1.0 / a, -(r_s / r).sqrt(), 0.0, 0.0);
        let e_r = FourVector::new_spherical(0.0, -1.0, 0.0, 0.0);
        let e_theta = FourVector::new_spherical(0.0, 0.0, -1.0, 0.0);
        let e_phi = FourVector::new_spherical(0.0, 0.0, 0.0, -1.0);

        let geometry = Schwarzschild::new(r_s, 1e-5);
        let vectors = vec![e_t, e_r, e_theta, e_phi];
        let orthonormal_vectors = gram_schmidt(&geometry, &point_in_spherical, &vectors);

        println!("{:?}", orthonormal_vectors);

        let norm = geometry.inner_product(
            &point_in_spherical,
            &orthonormal_vectors[0],
            &orthonormal_vectors[0],
        );
        assert_abs_diff_eq!(norm, 1.0);
        assert_eq!(orthonormal_vectors.len(), 4);
        for i in 1..4 {
            let norm = geometry.inner_product(
                &point_in_spherical,
                &orthonormal_vectors[i],
                &orthonormal_vectors[i],
            );
            assert_abs_diff_eq!(norm, -1.0);
        }

        for i in 0..4 {
            for j in 0..4 {
                if i != j {
                    assert_abs_diff_eq!(
                        geometry.inner_product(
                            &point_in_spherical,
                            &orthonormal_vectors[i],
                            &orthonormal_vectors[j]
                        ),
                        0.0
                    );
                }
            }
        }
    }

    #[test]
    fn test_gram_schmidt_schwarzschild() {
        let r_s = 2.0;

        let position = Point::new_cartesian(0.0, 3.0, 2.0, 1.0);
        let point_in_spherical = cartesian_to_spherical(&position);
        println!("{:?}", point_in_spherical);
        let r = point_in_spherical.vector[1];
        let theta = point_in_spherical.vector[2];
        let phi = point_in_spherical.vector[3];

        let a = 1.0 - r_s / r;

        let e_t = FourVector::new_spherical(1.0 / a, -(r_s / r).sqrt(), 0.0, 0.0);
        let e_1 = FourVector::new_spherical(0.0, -1.0, 0.0, 0.0);
        let e_2 = FourVector::new_spherical(0.0, 0.0, -1.0, 0.0);
        let e_3 = FourVector::new_spherical(0.0, 0.0, 0.0, -1.0);

        let geometry = Schwarzschild::new(r_s, 1e-5);
        let vectors = vec![e_t, e_1, e_2, e_3];
        let orthonormal_vectors = gram_schmidt(&geometry, &point_in_spherical, &vectors);

        assert_eq!(orthonormal_vectors.len(), 4);

        let e_t_expected = FourVector::new_spherical(1.0 / a, -(r_s / r).sqrt(), 0.0, 0.0);
        let e_1_expected = FourVector::new_spherical((r_s / r).sqrt() / a, -1.0, 0.0, 0.0);
        let e_2_expected = FourVector::new_spherical(0.0, 0.0, -1.0 / r, 0.0);
        let e_3_expected = FourVector::new_spherical(0.0, 0.0, 0.0, -1.0 / (r * theta.sin()));

        println!(
            "{:?}",
            [e_t_expected, e_1_expected, e_2_expected, e_3_expected]
        );
        println!("{:?}", orthonormal_vectors);

        let geometry = EuclideanSpace::new();

        assert_eq!(orthonormal_vectors[0], e_t_expected);
        assert_eq!(orthonormal_vectors[1], e_1_expected);
        assert_eq!(orthonormal_vectors[2], e_2_expected);
        assert_eq!(orthonormal_vectors[3], e_3_expected);
    }
}
