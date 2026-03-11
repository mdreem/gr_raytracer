use crate::rendering::black_body_radiation::get_cie_xyz_of_black_body_redshifted;
use crate::rendering::color::xyz_to_srgb;

pub fn run_blackbody(temperature: f64, redshift: f64) {
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

pub fn run_blackbody_spectrum(
    min_temperature: f64,
    max_temperature: f64,
    min_redshift: f64,
    max_redshift: f64,
    width: u32,
    height: u32,
    filename: String,
) {
    use indicatif::{ProgressBar, ProgressStyle};
    use rayon::prelude::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let total_pixels = (width * height) as usize;
    let mut pixels = vec![0u8; total_pixels * 4];

    let pb = ProgressBar::new(total_pixels as u64);
    pb.set_style(ProgressStyle::with_template("🔥 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
        .expect("Failed to create progress style")
        .progress_chars("█▇▆▅▄▃▂▁  "));

    let count = AtomicUsize::new(0);

    pixels.par_chunks_mut(4).enumerate().for_each(|(i, p)| {
        count.fetch_add(1, Ordering::SeqCst);
        pb.set_position(count.load(Ordering::Relaxed) as u64);

        let x = (i % width as usize) as f64;
        let y = (i / width as usize) as f64;

        let temperature = min_temperature
            + x * (max_temperature - min_temperature) / (width as f64 - 1.0);
        let redshift =
            min_redshift + y * (max_redshift - min_redshift) / (height as f64 - 1.0);

        let cie_xyz = get_cie_xyz_of_black_body_redshifted(temperature * redshift);
        let exposure = 1.0 / (cie_xyz.x + cie_xyz.y + cie_xyz.z);
        let color = xyz_to_srgb(&cie_xyz, exposure);

        p[0] = color.r;
        p[1] = color.g;
        p[2] = color.b;
        p[3] = color.alpha;
    });

    pb.finish();

    let imgbuf = image::RgbaImage::from_raw(width, height, pixels)
        .expect("Failed to create image buffer from raw pixels");
    imgbuf.save(&filename).expect("Failed to save spectrum image");
    println!("Saved blackbody spectrum to {}", filename);
}
