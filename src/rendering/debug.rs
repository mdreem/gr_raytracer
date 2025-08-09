use crate::geometry::geometry::Geometry;
use crate::geometry::point::Point;
use crate::rendering::camera::Camera;
use std::fs::File;
use std::io::Write;

pub fn save_rays_to_file<G: Geometry>(
    rows: i64,
    cols: i64,
    position: &Point,
    geometry: G,
    camera: Camera,
) {
    let mut file = File::create("rays.csv").expect("Unable to create file");
    file.write_all(b"row,col,m_x0,m_x1,m_x2,m_x3,m_sc\n")
        .expect("Unable to write file");

    for i in 1..rows {
        for j in 1..cols {
            let ray = camera.get_ray_for(i, j);
            let momentum = ray.momentum.get_as_vector();
            file.write_all(
                format!(
                    "{},{},{},{},{},{},{}\n",
                    i,
                    j,
                    momentum[0],
                    momentum[1],
                    momentum[2],
                    momentum[3],
                    geometry.inner_product(position, &ray.momentum, &ray.momentum),
                )
                .as_bytes(),
            )
            .expect("Unable to write file");
        }
    }
    println!("Finished writing rays.");
}
