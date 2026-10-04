use crate::geometry::coordinates::atlas::{Atlas, Chart};
use crate::geometry::point::Point;

#[derive(Copy, Clone, PartialEq, Debug)]
pub struct RotatedSpherical {
    pub rotation: f64,
}

pub struct SphericalAtlas {}

impl Atlas for RotatedSpherical {
    type Chart = SphericalChart;

    fn best_chart(&self, position: &Point<SphericalChart>, current: &Self::Chart) -> Self::Chart {
        let SphericalChart::Rotated(chart_data) = current;
        *current
    }

    fn transition_point(
        &self,
        position: &Point<SphericalChart>,
        from: &Self::Chart,
        to: &Self::Chart,
    ) -> Point<SphericalChart> {
        let SphericalChart::Rotated(from_rotation) = from;
        let SphericalChart::Rotated(to_rotation) = to;

        todo!()
    }

    fn transition_four_vector(
        &self,
        position: &Point<SphericalChart>,
        from: &Self::Chart,
        to: &Self::Chart,
    ) -> Point<SphericalChart> {
        todo!()
    }
}

impl Chart for SphericalChart {
    fn tetrad_handedness_fallback(&self) -> f64 {
        -1.0
    }
}

#[derive(Copy, Clone, PartialEq, Debug)]
pub enum SphericalChart {
    Rotated(RotatedSpherical),
}
