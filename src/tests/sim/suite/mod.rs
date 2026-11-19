//! Scenario based test suite.
//!
//! Each test describes a [Scenario], runs it through the public `Solver` API
//! and asserts on performance metrics against the simulated truth.
//! Thresholds are derived from the physics of the scenario (noise level,
//! geometry, unmodeled effects) with a safety margin, so that a regression
//! of the navigation filter or of a physical model is caught.
//!
//! Tests marked `#[ignore]` describe *targets* of upcoming milestones (they
//! fail today on purpose): run them with `cargo test -- --ignored`.
//!
//! Set `GNSS_RTK_SIM_OUTPUT=<dir>` to export every run as CSV.
mod absolute;
mod benchmark;
mod kinematic;
mod model;
mod ppp;
mod robustness;
mod rtk;

use crate::{
    prelude::{ClockProfile, Config, Frame, UserParameters, UserProfile},
    tests::{earth_frame, init_logger, sim::RunReport},
};

pub(crate) fn frame() -> Frame {
    init_logger();
    earth_frame()
}

pub(crate) fn params(profile: UserProfile) -> UserParameters {
    UserParameters::new(profile, ClockProfile::Quartz)
}

/// Disables the postfit smoother, to assess the navigation filter itself
pub(crate) fn no_postfit(cfg: &mut Config) {
    cfg.solver.postfit_denoising = 0.0;
}

/// Asserts with a readable report on failure
pub(crate) fn check(report: &RunReport, condition: bool, what: &str) {
    assert!(
        condition,
        "{} - {what}\n{report}\nfirst failures: {:?}",
        report.name,
        report.failures().iter().take(3).collect::<Vec<_>>()
    );
}
