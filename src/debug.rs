use crate::camera::Camera;
use crate::geometry::Geometry;
use nalgebra::Vector4;
use std::fs::File;
use std::io::Write;

pub fn save_rays_to_file<G: Geometry>(
    rows: i64,
    cols: i64,
    position: &Vector4<f64>,
    geometry: G,
    camera: Camera,
) {
    let mut file = File::create("rays.csv").expect("Unable to create file");
    file.write_all(b"row,col,d_x0,d_x1,d_x2,d_x3,m_x0,m_x1,m_x2,m_x3,m_sc,d_sc\n")
        .expect("Unable to write file");

    for i in 1..rows {
        for j in 1..cols {
            let ray = camera.get_ray_for(i, j);
            let dir = ray.direction.get_as_vector();
            let momentum = ray.momentum.get_as_vector();
            file.write_all(
                format!(
                    "{},{},{},{},{},{},{},{},{},{},{},{}\n",
                    i,
                    j,
                    dir[0],
                    dir[1],
                    dir[2],
                    dir[3],
                    momentum[0],
                    momentum[1],
                    momentum[2],
                    momentum[3],
                    geometry.mul(&position, &ray.momentum, &ray.momentum),
                    geometry.mul(&position, &ray.direction, &ray.direction),
                )
                .as_bytes(),
            )
            .expect("Unable to write file");
        }
    }
    println!("Finished writing rays.");
}
