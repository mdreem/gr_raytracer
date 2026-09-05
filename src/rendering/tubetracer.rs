use crate::rendering::scene::RaySample;

pub struct SampleTube<'a> {
    pub a: &'a RaySample,
    pub b: &'a RaySample,
    pub c: &'a RaySample,
    pub d: &'a RaySample,
}

impl SampleTube<'_> {
    pub fn new<'a>(
        a: &'a RaySample,
        b: &'a RaySample,
        c: &'a RaySample,
        d: &'a RaySample,
    ) -> SampleTube<'a> {
        SampleTube { a, b, c, d }
    }

    pub fn from_buffer(
        buffer: &[RaySample],
        row: u32,
        col: u32,
        width: u32,
    ) -> Option<SampleTube<'_>> {
        if col >= width - 1 {
            return None;
        }

        let idx_a = (row * width + col) as usize;
        let idx_b = (row * width + col + 1) as usize;
        let idx_c = ((row + 1) * width + col) as usize;
        let idx_d = ((row + 1) * width + col + 1) as usize;

        if !(idx_d < buffer.len()) {
            return None;
        }
        Some(SampleTube::new(
            &buffer[idx_a],
            &buffer[idx_b],
            &buffer[idx_c],
            &buffer[idx_d],
        ))
    }
}
