//! Precise Point Positioning targets.
//!
//! These scenarios describe the expected behavior of the PPP milestones
//! (float PPP: M2-M4, PPP-AR: M5). They are ignored until the corresponding
//! milestone lands and act as executable specifications.
use crate::{
    prelude::{Method, UserProfile},
    tests::sim::{IonoModel, NoiseModel, RunOptions, Scenario, Simulator, Technique, run},
};

use super::{check, frame, no_postfit, params};

fn static_ppp() -> RunOptions {
    RunOptions::new(Method::PPP, Technique::Absolute)
        .with_params(params(UserProfile::Static))
        .with_config(no_postfit)
}

/// Float PPP, noise free with ionosphere: mm level after initialization.
#[test]
#[ignore]
fn target_ideal_float_ppp() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::ideal_static("target_ideal_float_ppp")
            .with_iono(IonoModel::perfectly_modeled(20.0)),
        frame,
    );

    let report = run(&sim, &static_ppp());

    check(&report, report.availability() == 1.0, "availability");
    check(&report, report.stats().rms_3d_m < 0.01, "3D RMS");
}

/// Float PPP, geodetic receiver: classical PPP convergence (< 10 cm within 30').
#[test]
#[ignore]
fn target_realistic_float_ppp_convergence() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::realistic_static("target_realistic_float_ppp")
            .with_duration(3.0 * 3600.0, 30.0)
            .with_noise(NoiseModel::geodetic()),
        frame,
    );

    let report = run(&sim, &static_ppp());

    check(
        &report,
        report
            .convergence_time_s(0.1)
            .map(|t| t < 1800.0)
            .unwrap_or(false),
        "convergence below 10 cm within 30'",
    );
    check(
        &report,
        report.stats_from(3600.0).rms_3d_m < 0.05,
        "3D RMS after 1 h",
    );
    check(&report, report.clock_rms_ns_from(3600.0) < 0.3, "timing");
}
