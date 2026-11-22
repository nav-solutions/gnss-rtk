//! Precise Point Positioning (float PPP) scenarios.
//!
//! One navigation filter processes the iono free code and phase of every
//! satellite: the solution starts at the CPP level and converges to the PPP
//! accuracy as the geometry changes (no mode switch). These scenarios cover the
//! convergence, the phase tracking (rising / setting satellites, cycle slips,
//! data gaps), the troposphere and the receiver clock estimation.
//!
//! PPP-AR (integer ambiguities) remains a target (see `target_ppp_ar_*`).
use crate::{
    prelude::{Constellation, Method, SV, UserProfile},
    tests::{
        sim::{
            ClockModel, CycleSlip, IonoModel, Outage, RunOptions, RunReport, Scenario,
            SimConstellation, Simulator, Station, Technique, Trajectory, TropoModel, run,
        },
        spec::{Preset, assert_contract},
    },
};

use super::{check, frame, no_postfit, params};

/// Convergence period of the float PPP (geodetic receiver, 30 s sampling)
const CONVERGENCE_S: f64 = 1800.0;

fn static_ppp() -> RunOptions {
    RunOptions::new(Method::PPP, Technique::Absolute)
        .with_params(params(UserProfile::Static))
        .with_config(no_postfit)
}

fn static_cpp() -> RunOptions {
    RunOptions::new(Method::CPP, Technique::Absolute)
        .with_params(params(UserProfile::Static))
        .with_config(no_postfit)
}

/// Geodetic receiver, 3 h static, GPS + Galileo.
fn geodetic(name: &str) -> Scenario {
    Scenario::realistic_static(name)
        .with_duration(3.0 * 3600.0, 30.0)
        .with_constellations(&[SimConstellation::GPS, SimConstellation::Galileo])
}

/// Satellites in view at sample `k`, sorted by decreasing elevation
fn satellites_in_view(sim: &Simulator, k: usize) -> Vec<SV> {
    let mut svs = sim
        .observe(Station::Rover, k)
        .iter()
        .filter_map(|cd| {
            sim.signal_truth(Station::Rover, cd.sv, k)
                .map(|truth| (cd.sv, truth.elevation_deg))
        })
        .collect::<Vec<_>>();

    svs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    svs.into_iter().map(|(sv, _)| sv).collect()
}

fn summary(report: &RunReport) -> String {
    let s = report.stats_from(CONVERGENCE_S);
    format!(
        "{}: avail={:.3} conv(10cm)={:?} s, after {CONVERGENCE_S} s: rms H={:.3} U={:.3} 3D={:.3} max3D={:.3} m, clock rms={:.3} ns",
        report.name,
        report.availability(),
        report.convergence_time_s(0.1),
        s.rms_h_m,
        s.rms_u_m,
        s.rms_3d_m,
        s.max_3d_m,
        report.clock_rms_ns_from(CONVERGENCE_S),
    )
}

/// Float PPP, noise free with ionosphere: mm level from the first epochs.
#[test]
fn ppp_ideal_float() {
    let sim = Simulator::new(
        Scenario::ideal_static("ppp_ideal_float").with_iono(IonoModel::perfectly_modeled(20.0)),
        frame(),
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    check(&report, report.availability() == 1.0, "availability");
    check(&report, report.stats().rms_3d_m < 0.01, "3D RMS");
    check(&report, report.clock_rms_ns_from(0.0) < 0.05, "clock");
    assert_contract(Preset::PPP, &report, 0.0);
}

/// Float PPP, geodetic receiver, Galileo only: classical PPP convergence.
#[test]
fn ppp_realistic_galileo_convergence() {
    let sim = Simulator::new(
        Scenario::realistic_static("ppp_realistic_galileo").with_duration(3.0 * 3600.0, 30.0),
        frame(),
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    check(&report, report.availability() == 1.0, "availability");
    check(
        &report,
        report.convergence_time_s(0.1).is_some_and(|t| t < 2700.0),
        "convergence below 10 cm within 45'",
    );
    check(
        &report,
        report.stats_from(3600.0).rms_3d_m < 0.05,
        "3D RMS after 1 h",
    );
    check(&report, report.clock_rms_ns_from(3600.0) < 0.3, "timing");
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}

/// GPS + Galileo: twice as many satellites, faster convergence than Galileo only.
#[test]
fn ppp_multi_constellation_converges_faster() {
    let frame = frame();

    let galileo = run(
        &Simulator::new(
            Scenario::realistic_static("ppp_galileo").with_duration(3.0 * 3600.0, 30.0),
            frame,
        ),
        &static_ppp(),
    );

    let both = run(
        &Simulator::new(geodetic("ppp_gps_galileo"), frame),
        &static_ppp(),
    );

    println!("{}\n{}", summary(&galileo), summary(&both));

    let t_galileo = galileo.convergence_time_s(0.1).unwrap_or(f64::INFINITY);
    let t_both = both.convergence_time_s(0.1).unwrap_or(f64::INFINITY);

    check(
        &both,
        t_both < 1500.0,
        "GPS + Galileo convergence within 25'",
    );
    check(&both, t_both < t_galileo, "faster than Galileo only");
    check(
        &both,
        both.stats_from(3600.0).rms_3d_m < 0.03,
        "3D RMS after 1 h",
    );
    assert_contract(Preset::PPP, &both, CONVERGENCE_S);
}

/// After convergence, float PPP is (at least) an order of magnitude better than CPP.
#[test]
fn ppp_is_ten_times_better_than_cpp() {
    let sim = Simulator::new(geodetic("ppp_vs_cpp"), frame());

    let ppp = run(&sim, &static_ppp());
    let cpp = run(&sim, &static_cpp());

    let (ppp_rms, cpp_rms) = (
        ppp.stats_from(3600.0).rms_3d_m,
        cpp.stats_from(3600.0).rms_3d_m,
    );

    println!("{}\nCPP 3D RMS after 1 h: {cpp_rms:.3} m", summary(&ppp));

    check(&ppp, ppp_rms * 10.0 < cpp_rms, "PPP 10x better than CPP");

    // the very first solution is the code (CPP) initialization: same accuracy
    let first = |r: &RunReport| r.outcomes[0].error_3d_m().unwrap();
    check(
        &ppp,
        (first(&ppp) - first(&cpp)).abs() < 0.5,
        "PPP starts at the CPP accuracy",
    );
}

/// Long static session: satellites rise and set, phase arcs start and end,
/// the accuracy must not degrade.
#[test]
fn ppp_rising_and_setting_satellites() {
    let sim = Simulator::new(
        geodetic("ppp_rising_setting").with_duration(6.0 * 3600.0, 60.0),
        frame(),
    );

    // the session must actually see satellites rise and set
    let first = satellites_in_view(&sim, 0);
    let last = satellites_in_view(&sim, sim.scenario().num_epochs() - 1);
    let lost = first.iter().filter(|sv| !last.contains(sv)).count();
    let new = last.iter().filter(|sv| !first.contains(sv)).count();
    assert!(
        lost >= 4 && new >= 4,
        "{lost} set, {new} rose: invalid scenario"
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    check(&report, report.availability() == 1.0, "availability");

    let second_half = report.stats_from(3.0 * 3600.0);
    check(&report, second_half.rms_3d_m < 0.03, "3D RMS, last 3 h");
    check(&report, second_half.max_3d_m < 0.08, "max 3D, last 3 h");
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}

/// Cycle slips (on one carrier, both carriers, wide lane neutral, large) after
/// convergence: detected, the arcs restart, the solution is not degraded.
#[test]
fn ppp_cycle_slips() {
    let frame = frame();
    let clean = geodetic("ppp_no_slip");
    let sim = Simulator::new(clean.clone(), frame);

    let k = (5400.0 / clean.interval_s) as usize;
    let svs = satellites_in_view(&sim, k);
    assert!(svs.len() >= 8);

    let mut slipped = clean
        .clone()
        .with_cycle_slip(CycleSlip {
            sv: svs[0],
            at_s: 5400.0,
            cycles_1: 1,
            cycles_2: 0,
        })
        .with_cycle_slip(CycleSlip {
            sv: svs[1],
            at_s: 5400.0,
            cycles_1: 0,
            cycles_2: -1,
        })
        .with_cycle_slip(CycleSlip {
            sv: svs[2],
            at_s: 6000.0,
            cycles_1: 1,
            cycles_2: 1,
        })
        .with_cycle_slip(CycleSlip {
            sv: svs[3],
            at_s: 6600.0,
            cycles_1: 9,
            cycles_2: 7,
        })
        .with_cycle_slip(CycleSlip {
            sv: svs[4],
            at_s: 7200.0,
            cycles_1: 1000,
            cycles_2: -250,
        })
        .with_cycle_slip(CycleSlip {
            sv: svs[5],
            at_s: 7200.0,
            cycles_1: -5,
            cycles_2: 0,
        });

    slipped.name = "ppp_cycle_slips".to_string();

    let reference = run(&sim, &static_ppp());
    let report = run(&Simulator::new(slipped, frame), &static_ppp());

    println!("{}\n{}", summary(&reference), summary(&report));

    check(&report, report.availability() == 1.0, "availability");

    let after = report.stats_from(5400.0);
    let nominal = reference.stats_from(5400.0);

    check(
        &report,
        after.rms_3d_m < nominal.rms_3d_m + 0.01,
        "RMS not degraded by the slips",
    );
    check(&report, after.max_3d_m < 0.05, "no jump at the slips");
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}

/// Signal outages: short (arc continues) and long (arc restarts), then a
/// complete data gap: the filter carries the position across.
#[test]
fn ppp_outages_and_data_gap() {
    let frame = frame();
    let scenario = geodetic("ppp_outages");
    let sim = Simulator::new(scenario.clone(), frame);

    let k = (4000.0 / scenario.interval_s) as usize;
    let svs = satellites_in_view(&sim, k);

    let mut degraded = scenario
        // short masking: within the gap tolerance
        .with_outage(Outage {
            sv: svs[0],
            start_s: 4000.0,
            end_s: 4060.0,
        })
        // long masking: the arc restarts
        .with_outage(Outage {
            sv: svs[1],
            start_s: 4000.0,
            end_s: 4900.0,
        })
        .with_outage(Outage {
            sv: svs[2],
            start_s: 4500.0,
            end_s: 5400.0,
        });

    // complete data gap (receiver reset, tunnel...): 5'
    for sv in satellites_in_view(&sim, (6000.0 / 30.0) as usize) {
        degraded = degraded.with_outage(Outage {
            sv,
            start_s: 6000.0,
            end_s: 6300.0,
        });
    }

    let report = run(&Simulator::new(degraded, frame), &static_ppp());
    println!("{}", summary(&report));

    let gap = report
        .outcomes
        .iter()
        .filter(|o| o.t_s >= 6000.0 && o.t_s < 6300.0)
        .collect::<Vec<_>>();

    assert!(
        gap.iter().all(|o| o.solution.is_none()),
        "no solution expected without data"
    );

    // the position survives the gap (static): no new convergence needed
    check(
        &report,
        report.stats_from(6300.0).max_3d_m < 0.10,
        "accuracy after the data gap",
    );
    check(
        &report,
        report.stats_from(CONVERGENCE_S).rms_3d_m < 0.05,
        "3D RMS",
    );
}

/// Unmodeled zenith wet delay (10 cm): estimated, the height is not biased.
#[test]
fn ppp_estimates_the_wet_troposphere() {
    let sim = Simulator::new(
        geodetic("ppp_zwd").with_tropo(TropoModel {
            residual_zwd_m: 0.10,
        }),
        frame(),
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    let s = report.stats_from(3600.0);

    // an unestimated 10 cm ZWD maps into a ~20-30 cm height bias
    check(&report, s.mean_enu_m.2.abs() < 0.04, "no height bias");
    check(&report, s.rms_3d_m < 0.05, "3D RMS after 1 h");
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}

/// Free running receiver clock (large drift, random walk): estimated as a
/// white noise process every epoch, the timing follows at the phase level.
#[test]
fn ppp_free_running_clock() {
    let sim = Simulator::new(
        geodetic("ppp_free_clock").with_rover_clock(ClockModel {
            offset_s: 2.0E-4,
            drift_s_s: 5.0E-8,
            random_walk_s_sqrt_s: 2.0E-9,
        }),
        frame(),
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    check(&report, report.availability() == 1.0, "availability");
    check(&report, report.clock_rms_ns_from(3600.0) < 0.3, "timing");
    check(
        &report,
        report.stats_from(3600.0).rms_3d_m < 0.05,
        "3D RMS after 1 h",
    );
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}

/// Kinematic float PPP (car on a ring road, 20 m/s): the phase tracks the
/// trajectory, the accuracy after convergence stays at the decimeter level.
#[test]
fn ppp_kinematic() {
    let (lat, lon, h) = crate::tests::sim::scenario::REFERENCE_SITE_LAT_LON_H;
    let center = Trajectory::static_geodetic(lat, lon, h).initial_position_m();

    let sim = Simulator::new(
        geodetic("ppp_kinematic").with_trajectory(Trajectory::Circle {
            center_ecef_m: center,
            radius_m: 500.0,
            period_s: 157.0,
        }),
        frame(),
    );

    let opts = RunOptions::new(Method::PPP, Technique::Absolute)
        .with_params(params(UserProfile::Car))
        .with_config(no_postfit);

    let report = run(&sim, &opts);
    let cpp = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Car))
            .with_config(no_postfit),
    );

    println!(
        "{}\nCPP rms 3D after 1 h: {:.3} m",
        summary(&report),
        cpp.stats_from(3600.0).rms_3d_m
    );

    check(&report, report.availability() == 1.0, "availability");
    check(
        &report,
        report.stats_from(3600.0).rms_3d_m < 0.5 * cpp.stats_from(3600.0).rms_3d_m,
        "kinematic PPP better than kinematic CPP",
    );
    assert_contract(Preset::PPP, &report, 3600.0);
}

/// Target (PPP-AR): integer ambiguity resolution requires satellite phase
/// biases (not part of the [crate::prelude::SpacebornBias] interface yet).
/// Expected: cm level within 15'.
#[test]
#[ignore]
fn target_ppp_ar_fast_convergence() {
    let sim = Simulator::new(geodetic("target_ppp_ar"), frame());
    let report = run(&sim, &static_ppp());

    check(
        &report,
        report.convergence_time_s(0.02).is_some_and(|t| t < 900.0),
        "2 cm within 15'",
    );
}

/// GPS only (L1 / L2): same behavior as Galileo only (E1 / E5b).
#[test]
fn ppp_gps_only() {
    let sim = Simulator::new(
        geodetic("ppp_gps_only").with_constellations(&[SimConstellation::GPS]),
        frame(),
    );

    let report = run(&sim, &static_ppp());
    println!("{}", summary(&report));

    // GPS only (24 satellites): a few epochs exceed the GDOP limit
    check(&report, report.availability() > 0.95, "availability");
    check(
        &report,
        report.stats_from(3600.0).rms_3d_m < 0.05,
        "3D RMS after 1 h",
    );
    assert!(
        (0..sim.scenario().num_epochs()).all(|k| sim
            .observe(Station::Rover, k)
            .iter()
            .all(|cd| cd.sv.constellation == Constellation::GPS)),
        "GPS only scenario"
    );
    assert_contract(Preset::PPP, &report, CONVERGENCE_S);
}
