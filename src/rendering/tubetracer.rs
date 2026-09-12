use crate::rendering::scene::RaySample;

/// A single ray tube's rectangular footprint on the camera screen.
///
/// Coordinates are absolute, continuous pixel-center coordinates: integer
/// rows/columns locate pixel centers, as in Camera::get_ray_for. These are
/// neither whole-screen extents nor world-space or distant-sky bounds.
/// Subdivision partitions this rectangle; the traced corner rays determine
/// the tube's distant-sky footprint separately.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TubeScreenBounds {
    /// Minimum row coordinate.
    pub top: f64,
    /// Maximum row coordinate.
    pub bottom: f64,
    /// Minimum column coordinate.
    pub left: f64,
    /// Maximum column coordinate.
    pub right: f64,
}

/// Four traced corner samples and their originating camera-screen rectangle.
/// Corners a/b/c/d are top-left/top-right/bottom-left/bottom-right.
pub struct SampleTube<'a> {
    pub a: &'a RaySample,
    pub b: &'a RaySample,
    pub c: &'a RaySample,
    pub d: &'a RaySample,
    /// Actual sampling positions; Ray.row/col retain integer pixel identities.
    pub screen_bounds: TubeScreenBounds,
}

impl SampleTube<'_> {
    pub fn new<'a>(
        a: &'a RaySample,
        b: &'a RaySample,
        c: &'a RaySample,
        d: &'a RaySample,
        screen_bounds: TubeScreenBounds,
    ) -> SampleTube<'a> {
        SampleTube {
            a,
            b,
            c,
            d,
            screen_bounds,
        }
    }

    /// Center, top, left, right, bottom: each is traced once and shared by
    /// the children. None means the bounds cannot be refined in f64.
    pub fn subdivision_points(&self) -> Option<[(f64, f64); 5]> {
        let TubeScreenBounds {
            top,
            bottom,
            left,
            right,
        } = self.screen_bounds;
        let row = top + (bottom - top) * 0.5;
        let col = left + (right - left) * 0.5;
        let finite = [top, bottom, left, right, row, col]
            .iter()
            .all(|v| v.is_finite());
        if !(finite && top < row && row < bottom && left < col && col < right) {
            return None;
        }
        Some([
            (row, col),
            (top, col),
            (row, left),
            (row, right),
            (bottom, col),
        ])
    }

    /// The samples must follow subdivision_points' ordering. Borrowing the
    /// samples keeps shared child boundaries identical without retracing.
    /// Each child's screen bounds cover one quarter of the parent's rectangle.
    pub fn children<'b>(&'b self, samples: &'b [RaySample; 5]) -> [SampleTube<'b>; 4] {
        let [center, top_mid, left_mid, right_mid, bottom_mid] = samples;
        let TubeScreenBounds {
            top,
            bottom,
            left,
            right,
        } = self.screen_bounds;
        let row = top + (bottom - top) * 0.5;
        let col = left + (right - left) * 0.5;
        [
            SampleTube::new(
                self.a,
                top_mid,
                left_mid,
                center,
                TubeScreenBounds {
                    top,
                    bottom: row,
                    left,
                    right: col,
                },
            ),
            SampleTube::new(
                top_mid,
                self.b,
                center,
                right_mid,
                TubeScreenBounds {
                    top,
                    bottom: row,
                    left: col,
                    right,
                },
            ),
            SampleTube::new(
                left_mid,
                center,
                self.c,
                bottom_mid,
                TubeScreenBounds {
                    top: row,
                    bottom,
                    left,
                    right: col,
                },
            ),
            SampleTube::new(
                center,
                right_mid,
                bottom_mid,
                self.d,
                TubeScreenBounds {
                    top: row,
                    bottom,
                    left: col,
                    right,
                },
            ),
        ]
    }

    /// Builds a base tube spanning four neighboring pixel centers, not the
    /// edges of a single output pixel.
    pub fn from_buffer(
        buffer: &[RaySample],
        row: u32,
        col: u32,
        width: u32,
    ) -> Option<SampleTube<'_>> {
        if width < 2 || col >= width - 1 {
            return None;
        }

        let idx_a = (row * width + col) as usize;
        let idx_b = (row * width + col + 1) as usize;
        let idx_c = ((row + 1) * width + col) as usize;
        let idx_d = ((row + 1) * width + col + 1) as usize;

        if idx_d >= buffer.len() {
            return None;
        }
        Some(SampleTube::new(
            &buffer[idx_a],
            &buffer[idx_b],
            &buffer[idx_c],
            &buffer[idx_d],
            TubeScreenBounds {
                top: buffer[idx_a].ray.row as f64,
                bottom: buffer[idx_c].ray.row as f64,
                left: buffer[idx_a].ray.col as f64,
                right: buffer[idx_b].ray.col as f64,
            },
        ))
    }
}
