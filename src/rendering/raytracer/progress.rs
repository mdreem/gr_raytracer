//! Terminal and headless render progress reporting.

use crate::rendering::error::RaytracerError;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use std::io::{self, IsTerminal};
use std::sync::atomic::{AtomicU64, Ordering};

/// Render progress that adapts to where it runs: the indicatif bar on an
/// interactive terminal, and plain `progress: <label> NN% (done/total)` lines
/// when stderr is not a TTY (e.g. a headless render pod, where indicatif draws
/// nothing at all). The plain lines use eprintln, not the log crate, so they
/// are not gated by RUST_LOG.
pub(super) struct Progress {
    bar: ProgressBar,
    is_tty: bool,
    total: u64,
    label: &'static str,
    last_percent: AtomicU64,
}

impl Progress {
    pub(super) fn new(total: u64, label: &'static str) -> Result<Self, RaytracerError> {
        let is_tty = io::stderr().is_terminal();
        let bar = ProgressBar::new(total);
        if is_tty {
            bar.set_style(
                ProgressStyle::with_template("🎨 {spinner:.green} [{elapsed_precise}] [{wide_bar:.blue}] {pos}/{len} ({percent_precise}%, {eta})")
                    .map_err(RaytracerError::ProgressBarTemplateError)?
                    .progress_chars("█▇▆▅▄▃▂▁  "),
            );
        } else {
            bar.set_draw_target(ProgressDrawTarget::hidden());
        }
        Ok(Self {
            bar,
            is_tty,
            total,
            label,
            last_percent: AtomicU64::new(0),
        })
    }

    /// Record that `done` units are complete. Safe to call from many threads;
    /// in the plain path each whole percent is printed at most once.
    pub(super) fn set(&self, done: u64) {
        if self.is_tty {
            self.bar.set_position(done);
            return;
        }
        if self.total == 0 {
            return;
        }
        let percent = done * 100 / self.total;
        let prev = self.last_percent.load(Ordering::Relaxed);
        if percent > prev
            && self
                .last_percent
                .compare_exchange(prev, percent, Ordering::SeqCst, Ordering::Relaxed)
                .is_ok()
        {
            eprintln!(
                "progress: {} {}% ({}/{})",
                self.label, percent, done, self.total
            );
        }
    }

    pub(super) fn finish(&self) {
        if self.is_tty {
            self.bar.finish();
        } else {
            eprintln!(
                "progress: {} 100% ({}/{})",
                self.label, self.total, self.total
            );
        }
    }
}
