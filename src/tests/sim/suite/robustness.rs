//! Robustness: outages, clocks, receiver grades, cycle slips, long runs.
use crate::{
    prelude::{Constellation, Method, SV, UserProfile},
    tests::sim::{
        ClockModel, CycleSlip, NoiseModel, Outage, RunOptions, Scenario, Simulator, Technique, run,
    },
};

use crate::tests::spec::{Preset, assert_contract};

use super::{check, frame, no_postfit, params};

fn static_cpp() -> RunOptions {
    RunOptions::new(Method::CPP, Technique::Absolute)
        .with_params(params(UserProfile::Static))
        .with_config(no_postfit)
}

/// All signals lost for 5 minutes: no solution during the outage, the filter
/// must resume afterwards without degradation.
#[test]
fn total_signal_outage() {
    let frame = frame();

    let mut scenario = Scenario::realistic_static("total_outage");

    for prn in 1..=24 {
        scenario = scenario.with_outage(Outage {
            sv: SV::new(Constellation::Galileo, prn),
            start_s: 1200.0,
            end_s: 1500.0,
        });
    }

    let sim = Simulator::new(scenario, frame);
    let report = run(&sim, &static_cpp());

    let during = report
        .outcomes
        .iter()
        .filter(|o| o.t_s >= 1200.0 && o.t_s < 1500.0);

    for o in during {
        assert!(o.solution.is_none(), "no solution expected at t={}s", o.t_s);
    }

    let after = report
        .outcomes
        .iter()
        .filter(|o| o.t_s >= 1500.0)
        .all(|o| o.solution.is_some());

    check(&report, after, "solutions must resume after the outage");
    check(
        &report,
        report.stats_from(1800.0).rms_3d_m < 0.8,
        "accuracy after the outage",
    );
    assert_contract(Preset::CPP, &report, 1800.0);
}

/// Partial outage (2 satellites masked for 10 minutes, one high, one low).
#[test]
fn partial_outage() {
    let frame = frame();

    let mut scenario = Scenario::realistic_static("partial_outage");

    for prn in [19, 17] {
        scenario = scenario.with_outage(Outage {
            sv: SV::new(Constellation::Galileo, prn),
            start_s: 900.0,
            end_s: 1500.0,
        });
    }

    let sim = Simulator::new(scenario, frame);
    let report = run(&sim, &static_cpp());

    check(&report, report.availability() == 1.0, "availability");
    check(&report, report.stats_from(600.0).rms_3d_m < 0.8, "3D RMS");
    assert_contract(Preset::CPP, &report, 600.0);
}

/// Free running oscillator: large offset, large drift, random walk.
#[test]
fn free_running_receiver_clock() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::realistic_static("free_running_clock").with_rover_clock(ClockModel {
            offset_s: 3.0E-4,
            drift_s_s: 5.0E-8,
            random_walk_s_sqrt_s: 1.0E-8,
        }),
        frame,
    );

    let report = run(&sim, &static_cpp());

    check(&report, report.availability() == 1.0, "availability");
    check(&report, report.stats_from(600.0).rms_3d_m < 0.8, "3D RMS");
    check(
        &report,
        report.clock_rms_ns_from(600.0) < 4.0,
        "clock tracking",
    );
    assert_contract(Preset::CPP, &report, 600.0);
}

/// Low cost receiver, single frequency.
#[test]
fn low_cost_receiver_spp() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::realistic_static("low_cost_spp").with_noise(NoiseModel::low_cost()),
        frame,
    );

    let report = run(
        &sim,
        &RunOptions::new(Method::SPP, Technique::Absolute)
            .with_params(params(UserProfile::Static))
            .with_config(no_postfit),
    );

    check(&report, report.availability() == 1.0, "availability");
    // dominated by the 50 % ionosphere model error (~3 m vertical)
    check(&report, report.stats_from(600.0).rms_3d_m < 5.0, "3D RMS");
    check(
        &report,
        report.stats_from(600.0).rms_h_m < 1.0,
        "horizontal RMS",
    );
}

/// Code based methods must not be affected by cycle slips.
#[test]
fn cycle_slips_do_not_affect_code_navigation() {
    let frame = frame();

    let clean = Scenario::realistic_static("no_slips").with_duration(1800.0, 30.0);

    let slipped = clean
        .clone()
        .with_cycle_slip(CycleSlip {
            sv: SV::new(Constellation::Galileo, 2),
            at_s: 600.0,
            cycles_1: 1000,
            cycles_2: -250,
        })
        .with_cycle_slip(CycleSlip {
            sv: SV::new(Constellation::Galileo, 18),
            at_s: 900.0,
            cycles_1: 1,
            cycles_2: 1,
        });

    let a = run(&Simulator::new(clean, frame), &static_cpp());
    let b = run(&Simulator::new(slipped, frame), &static_cpp());

    assert_eq!(a.stats().rms_3d_m, b.stats().rms_3d_m);
}

/// Six hours: satellites rise and set, the solution must stay available and stable.
#[test]
fn six_hours_static() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::realistic_static("six_hours").with_duration(6.0 * 3600.0, 60.0),
        frame,
    );

    let report = run(&sim, &static_cpp());
    let stats = report.stats_from(1800.0);

    check(&report, report.availability() == 1.0, "availability");
    check(&report, stats.rms_3d_m < 0.5, "3D RMS");
    check(&report, stats.max_3d_m < 1.5, "3D max");
    assert_contract(Preset::CPP, &report, 1800.0);
}
