//! Differential navigation (`Solver::rtk`).
use nalgebra::Vector3;

use crate::{
    prelude::{Method, Rc, Solver, UserProfile},
    tests::{
        almanac,
        ephemeris::NullEph,
        sim::{IonoModel, RunOptions, Scenario, Simulator, Technique, run},
        time::NullTime,
    },
};

use super::{check, frame, no_postfit, params};

/// Noise free RTK must be exact for short baselines. This used to fail by
/// meters (far field linearization of the double differences).
#[test]
fn ideal_rtk_short_baselines() {
    let frame = frame();

    for (baseline_m, max_rms_m) in [(1_000.0, 0.005), (10_000.0, 0.04)] {
        let sim = Simulator::new(
            Scenario::ideal_static(&format!("ideal_rtk_{baseline_m}")).with_base_offset_enu(
                0.6 * baseline_m,
                0.8 * baseline_m,
                0.0,
            ),
            frame,
        );

        for method in [Method::SPP, Method::CPP] {
            let report = run(&sim, &RunOptions::new(method, Technique::Rtk));
            check(&report, report.availability() == 1.0, "availability");
            check(&report, report.stats().rms_3d_m < max_rms_m, "3D RMS");
        }
    }
}

/// Long baselines: the (unmodeled) differential troposphere grows with the
/// baseline, but the solution must remain decimetric up to 30 km.
#[test]
fn ideal_rtk_long_baseline() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::ideal_static("ideal_rtk_30km").with_base_offset_enu(18_000.0, 24_000.0, 0.0),
        frame,
    );

    let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Rtk));
    check(&report, report.stats().rms_3d_m < 0.2, "30 km baseline");
}

/// Short baseline: the ionosphere cancels in the double differences, even
/// in single frequency and without any ionosphere model.
#[test]
fn rtk_spp_cancels_the_ionosphere() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::ideal_static("rtk_spp_iono")
            .with_iono(IonoModel {
                vtec_tecu: 50.0,
                model_error: 1.0,
            })
            .with_base_offset_enu(1_500.0, 2_000.0, 0.0),
        frame,
    );

    let report = run(&sim, &RunOptions::new(Method::SPP, Technique::Rtk));
    check(&report, report.stats().rms_3d_m < 0.05, "iono must cancel");
}

/// Realistic noise, static rover, 5 km baseline.
#[test]
fn realistic_static_rtk() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::realistic_static("realistic_static_rtk")
            .with_base_offset_enu(3_000.0, 4_000.0, 0.0),
        frame,
    );

    for (method, max_rms_m) in [(Method::SPP, 0.35), (Method::CPP, 1.0)] {
        let report = run(
            &sim,
            &RunOptions::new(method, Technique::Rtk)
                .with_params(params(UserProfile::Static))
                .with_config(no_postfit),
        );

        check(&report, report.availability() == 1.0, "availability");
        check(
            &report,
            report.stats_from(600.0).rms_3d_m < max_rms_m,
            "3D RMS after 10'",
        );
    }
}

/// The README promises that one may switch between RTK and PPP at any epoch
/// (e.g. RTK network loss), the filter continuing from its past state.
#[test]
fn context_switch_rtk_ppp_rtk() {
    let frame = frame();

    let scenario = Scenario::ideal_static("context_switch")
        .with_duration(1800.0, 30.0)
        .with_base_offset_enu(3_000.0, 4_000.0, 0.0);

    let sim = Simulator::new(scenario.clone(), frame);

    let opts = RunOptions::new(Method::CPP, Technique::Rtk);

    let mut solver = Solver::new(
        almanac(),
        frame,
        opts.cfg.clone(),
        Rc::new(NullEph {}),
        Rc::new(sim.orbit_source()),
        Rc::new(sim.spaceborn_biases()),
        Rc::new(sim.environment()),
        NullTime {},
        None,
    );

    let base = sim.rtk_base();

    for epoch in sim.rover_epochs() {
        // RTK network lost between 10' and 20'
        let network_lost = epoch.t_s >= 600.0 && epoch.t_s < 1200.0;

        let pvt = if network_lost {
            solver.ppp(epoch.epoch, opts.params, &epoch.candidates)
        } else {
            solver.rtk(epoch.epoch, opts.params, &epoch.candidates, &base)
        }
        .unwrap_or_else(|e| panic!("t={}s (network lost={network_lost}): {e}", epoch.t_s));

        let est = Vector3::new(pvt.pos_m.0, pvt.pos_m.1, pvt.pos_m.2);
        let err = (est - epoch.truth.position_m).norm();

        assert!(
            err < 0.05,
            "t={}s (network lost={network_lost}): error {err} m",
            epoch.t_s
        );
    }
}
