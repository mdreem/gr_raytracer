use crate::cli::cli::GlobalOpts;
use crate::configuration::RenderConfig;
use crate::geometry::atlases::Atlas;
use crate::geometry::circular_orbit::OrbitKillingDecomposition;
use crate::geometry::coordinates::cartesian::CartesianChart;
use crate::geometry::coordinates::chart::{Chart, GeometryChart};
use crate::geometry::four_vector::FourVector;
use crate::geometry::point::Point;
use crate::geometry::tetrad::Tetrad;
use crate::rendering::error::RaytracerError;
use crate::rendering::integration::runge_kutta::OdeFunction;
use crate::rendering::integration::state::EquationOfMotionState;
use crate::rendering::light::temperature::TemperatureComputer;
use nalgebra::{Const, Matrix4};
use std::io::Write;

pub trait GeodesicSolver: OdeFunction<Const<8>> {
    type Chart: Chart;
    fn chart(&self) -> Self::Chart;
    fn geodesic(&self, t: f64, y: &EquationOfMotionState) -> EquationOfMotionState;
    fn create_initial_state(&self, momentum: &FourVector<Self::Chart>) -> EquationOfMotionState {
        EquationOfMotionState::from_column_slice(&[
            momentum.at[0],
            momentum.at[1],
            momentum.at[2],
            momentum.at[3],
            momentum[0],
            momentum[1],
            momentum[2],
            momentum[3],
        ])
    }
    fn momentum_from_state(&self, y: &EquationOfMotionState) -> FourVector<Self::Chart> {
        FourVector::new(
            y[4],
            y[5],
            y[6],
            y[7],
            Point::new(y[0], y[1], y[2], y[3], self.chart()),
        )
    }
}

pub trait InnerProduct: HasChart {
    fn inner_product(
        &self,
        position: &Point<Self::Chart>,
        v: &FourVector<Self::Chart>,
        w: &FourVector<Self::Chart>,
    ) -> f64;
}

pub trait Signature {
    fn signature(&self) -> [f64; 4];
}

/// Support quantities that can be derived from the geometry.
/// This includes things like computing various velocities at a given position or ways to
/// compute temperatures.
pub trait SupportQuantities: HasChart {
    /// The static (Killing) observer: at rest relative to infinity, purely
    /// along d_t. Does not exist inside the ergosphere of a spinning hole
    /// (components become non-finite there).
    fn get_stationary_velocity_at(&self, position: &Point<Self::Chart>) -> FourVector<Self::Chart>;
    /// The ZAMO (zero angular momentum observer / locally non-rotating
    /// frame): co-rotates with frame dragging, u . d_phi = 0. Coincides with
    /// the static observer wherever there is no rotation; exists everywhere
    /// outside the horizon, including inside the ergosphere.
    fn get_zamo_velocity_at(&self, position: &Point<Self::Chart>) -> FourVector<Self::Chart>;
    fn get_circular_orbit_velocity_at(
        &self,
        position: &Point<Self::Chart>,
    ) -> Result<FourVector<Self::Chart>, RaytracerError>;
    /// The axial Killing vector d_phi (generator of rotations about the spin
    /// axis) expressed in this geometry's chart: (0,0,0,1) in spherical/BL
    /// charts, (0, -y, x, 0) in Cartesian charts. Contracting the ray momentum
    /// with it yields the conserved p_phi.
    fn axial_killing_vector(&self, position: &Point<Self::Chart>) -> FourVector<Self::Chart>;
    /// Killing coefficients (u^t, u^phi) of the equatorial circular orbit at
    /// the position's radius. Chart-independent scalars; see
    /// `OrbitKillingDecomposition` for how they pair with conserved (p_t, p_phi).
    fn circular_orbit_killing_coefficients(
        &self,
        position: &Point<Self::Chart>,
    ) -> Result<OrbitKillingDecomposition, RaytracerError>;
    fn get_temperature_computer(
        &self,
        temperature: f64,
        inner_radius: f64,
        outer_radius: f64,
    ) -> Result<Box<dyn TemperatureComputer>, RaytracerError>;
}

#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ConstantsOfMotion {
    constants: Vec<(&'static str, f64)>,
}

#[allow(dead_code)]
impl ConstantsOfMotion {
    pub fn push(&mut self, name: &'static str, value: f64) {
        self.constants.push((name, value));
    }

    pub fn as_slice(&self) -> &[(&'static str, f64)] {
        &self.constants
    }
}

/// Radius, in units of the Schwarzschild radius `r_s`, within which a photon
/// that has exhausted the integrator's step budget without escaping or being
/// absorbed is classified as trapped (a bound / closed orbit). The threshold
/// comfortably contains the photon sphere (1.5·r_s) and the accretion disc
/// while staying well inside the celestial sphere. Shared by the
/// `closed_orbit` implementations of the black-hole geometries.
pub const TRAPPED_ORBIT_RADIUS_FACTOR: f64 = 5.0;

pub trait HasChart {
    type Chart: GeometryChart;
}

pub trait Geometry: HasChart + InnerProduct + SupportQuantities + Signature + Sync {
    fn get_tetrad_at(&self, position: &Point<Self::Chart>) -> Tetrad<Self::Chart>;
    fn lorentz_transformation(
        &self,
        position: &Point<Self::Chart>,
        velocity: &FourVector<Self::Chart>,
    ) -> Matrix4<f64>;
    fn inside_horizon(&self, position: &Point<Self::Chart>) -> bool;
    fn closed_orbit(
        &self,
        position: &Point<Self::Chart>,
        step_index: usize,
        max_steps: usize,
    ) -> bool;
    fn get_geodesic_solver(
        &self,
        momentum: &FourVector<Self::Chart>,
    ) -> Box<dyn GeodesicSolver<Chart = Self::Chart>>;
    fn get_radial_coordinate(&self, position: &Point<Self::Chart>) -> f64;
    /// The chart to continue integrating in at `position`: the atlas' choice by default.
    /// Geometries whose solver handles only some of their charts restrict it.
    fn integration_chart(&self, position: &Point<Self::Chart>) -> Self::Chart {
        position.chart.best_chart(position)
    }
    /// d(affine parameter)/d(`Step::t`) at `position`: 1 for solvers that integrate in the
    /// affine parameter, Σ for ones that integrate in Mino time (KerrBL).
    fn step_parameter_scale(&self, _position: &Point<Self::Chart>) -> f64 {
        1.0
    }
    /// Return the invariants along a geodesic, e.g. E and L_z.
    #[allow(dead_code)]
    fn get_constants_of_motion(
        &self,
        position: &Point<Self::Chart>,
        momentum: &FourVector<Self::Chart>,
    ) -> ConstantsOfMotion;
}

pub trait RenderableGeometry {
    fn render(
        &self,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point<CartesianChart>,
        filename: String,
        from_row: Option<u32>,
        from_col: Option<u32>,
        to_row: Option<u32>,
        to_col: Option<u32>,
    ) -> Result<crate::cli::shared::RenderReport, RaytracerError>;

    fn render_ray(
        &self,
        row: i64,
        col: i64,
        opts: GlobalOpts,
        config: RenderConfig,
        camera_position: Point<CartesianChart>,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError>;

    fn render_ray_at(
        &self,
        position: Point<CartesianChart>,
        direction: FourVector<CartesianChart>,
        opts: GlobalOpts,
        write: &mut dyn Write,
    ) -> Result<(), RaytracerError>;
}
