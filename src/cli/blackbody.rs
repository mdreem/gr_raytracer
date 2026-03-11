use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::xyz_to_srgb;

pub fn print_blackbody_color(temperature: f64, redshift: f64) {
    let cie_xyz = get_cie_xyz_of_black_body_redshifted(temperature * redshift);
    let exposure = 1.0 / (cie_xyz.x + cie_xyz.y + cie_xyz.z);
    let color = xyz_to_srgb(&cie_xyz, exposure);
    println!(
        "Blackbody color at T={}K (redshift={}):",
        temperature, redshift
    );
    println!("XYZ:  {:.4}, {:.4}, {:.4}", cie_xyz.x, cie_xyz.y, cie_xyz.z);
    println!("sRGB: R={}, G={}, B={}", color.r, color.g, color.b);
}

