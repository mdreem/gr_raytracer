use crate::geometry::coordinates::atlas::Chart;

pub struct KerrAtlas {}

/// Boyer-Lindquist coordinates (t, r, θ, φ) for the Kerr metric.
///
/// The Cartesian embedding used is the null-tetrad convention consistent with
/// the Kerr-Schild form used in `geometry/kerr.rs`:
///   x = (r cosφ − a sinφ) sinθ
///   y = (r sinφ + a cosφ) sinθ
///   z = r cosθ
///
/// Note: this differs from the oblate-spheroidal convention
/// (x = √(r²+a²) sinθ cosφ) used in some references.
///
/// # Note on equality
/// `PartialEq` is derived, which compares `a` and `r_s` with exact `f64` equality.
/// This is correct as long as they are always initialised from the same
/// values (e.g. `kerr_bl.a`) rather than an independent floating-point
/// computation.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BoyerLindquistData {
    pub a: f64,
    pub r_s: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct KerrSchildData {
    pub a: f64,
    pub r_s: f64,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum KerrChart {
    BoyerLindquist(BoyerLindquistData),
    KerrSchild(KerrSchildData),
}

impl Chart for KerrChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        match self {
            KerrChart::BoyerLindquist(_) => -1.0,
            KerrChart::KerrSchild(_) => 1.0,
        }
    }
}
