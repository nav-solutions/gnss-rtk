//! Absolute navigation (`Solver::ppp`) with code based methods (SPP, CPP).
use crate::{
    prelude::{Method, UserProfile},
    tests::sim::{
        ClockModel, Initialization, IonoModel, RunOptions, Scenario, SimConstellation, Simulator,
        Technique, TropoModel, run,
    },
};

use crate::tests::spec::{Preset, assert_contract};

use super::{check, frame, no_postfit, params};

/// Noise free, perfectly modeled: the whole chain must be exact (mm level).
#[test]
fn ideal_static_spp_cpp_survey() {
    let frame = frame();
    let sim = Simulator::new(Scenario::ideal_static("ideal_static"), frame);

    for method in [Method::SPP, Method::CPP] {
        let report = run(&sim, &RunOptions::new(method, Technique::Absolute));
        let stats = report.stats();

        check(&report, report.availability() == 1.0, "availability");
        check(&report, stats.rms_3d_m < 0.01, "3D RMS must be < 1 cm");
        check(&report, stats.max_3d_m < 0.02, "3D max must be < 2 cm");
        check(
            &report,
            report.clock_rms_ns_from(0.0) < 0.05,
            "clock RMS must be < 0.05 ns",
        );
        assert_contract(Preset::new(method, false), &report, 0.0);
    }
}

#[test]
fn ideal_static_multi_constellation() {
    let frame = frame();
    let sim = Simulator::new(
        Scenario::ideal_static("ideal_gps_gal")
            .with_constellations(&[SimConstellation::GPS, SimConstellation::Galileo]),
        frame,
    );

    let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Absolute));

    check(&report, report.availability() == 1.0, "availability");
    check(
        &report,
        report.stats().rms_3d_m < 0.02,
        "3D RMS must be < 2 cm",
    );
    assert_contract(Preset::CPP, &report, 0.0);
}

/// Wrong apriori (1 km): the solver must converge to the truth.
#[test]
fn ideal_static_wrong_apriori() {
    let frame = frame();
    let sim = Simulator::new(Scenario::ideal_static("wrong_apriori"), frame);

    let opts = RunOptions::new(Method::CPP, Technique::Absolute)
        .with_initialization(Initialization::Apriori(800.0, -500.0, 300.0))
        .with_params(params(UserProfile::Airplane))
        .with_config(no_postfit);

    let report = run(&sim, &opts);

    check(&report, report.availability() == 1.0, "availability");
    // first epoch suffers from models evaluated at the wrong apriori (troposphere)
    check(
        &report,
        report.stats_from(30.0).max_3d_m < 0.02,
        "must converge after 1 epoch",
    );
    assert_contract(Preset::CPP, &report, 30.0);
}

/// Regression test (Sagnac): a large receiver clock offset must not bias the solution
/// once the clock is known (geometric time of flight).
#[test]
fn receiver_clock_offset_does_not_bias() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::ideal_static("rx_clock_offset")
            .with_duration(600.0, 30.0)
            .with_rover_clock(ClockModel {
                offset_s: 5.0E-4,
                drift_s_s: 1.0E-9,
                random_walk_s_sqrt_s: 0.0,
            }),
        frame,
    );

    let opts = RunOptions::new(Method::CPP, Technique::Absolute)
        .with_params(params(UserProfile::Airplane))
        .with_config(no_postfit);

    let report = run(&sim, &opts);

    check(
        &report,
        report.stats_from(30.0).max_3d_m < 0.02,
        "Sagnac with receiver clock",
    );
    check(
        &report,
        report.clock_rms_ns_from(30.0) < 0.1,
        "clock estimate",
    );
    assert_contract(Preset::CPP, &report, 30.0);
}

/// Static, geodetic grade noise: CPP accumulates information with the static profile.
#[test]
fn realistic_static_cpp_convergence() {
    let frame = frame();
    let sim = Simulator::new(Scenario::realistic_static("realistic_static_cpp"), frame);

    let opts = RunOptions::new(Method::CPP, Technique::Absolute)
        .with_params(params(UserProfile::Static))
        .with_config(no_postfit);

    let report = run(&sim, &opts);
    let stats = report.stats_from(600.0);

    check(&report, report.availability() == 1.0, "availability");
    check(
        &report,
        stats.rms_3d_m < 0.8,
        "3D RMS after 10' must be < 0.8 m",
    );
    check(
        &report,
        stats.max_3d_m < 2.0,
        "3D max after 10' must be < 2 m",
    );
    check(
        &report,
        report
            .convergence_time_s(1.0)
            .map(|t| t < 1800.0)
            .unwrap_or(false),
        "must converge under 1 m within 30'",
    );
    check(
        &report,
        report.clock_rms_ns_from(600.0) < 4.0,
        "clock RMS < 4 ns",
    );
    assert_contract(Preset::CPP, &report, 600.0);
}

/// Regression test (Q from the user profile): accumulating the measurements of a
/// static receiver must beat an epoch-wise solution.
#[test]
fn static_profile_beats_pedestrian_profile() {
    let frame = frame();
    let sim = Simulator::new(Scenario::realistic_static("profile"), frame);

    let static_run = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Static))
            .with_config(no_postfit),
    );

    let pedestrian_run = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Pedestrian))
            .with_config(no_postfit),
    );

    let (s, p) = (
        static_run.stats_from(600.0),
        pedestrian_run.stats_from(600.0),
    );

    check(
        &static_run,
        s.rms_3d_m < 0.5 * p.rms_3d_m,
        "static must be at least 2x better",
    );
}

/// Regression test: the ionospheric model must not be applied to the iono free
/// combination. CPP results must not depend on the ionosphere model quality.
#[test]
fn cpp_is_independent_of_the_ionosphere_model() {
    let frame = frame();

    let mut rms = Vec::new();

    for (name, iono) in [
        ("iono_perfect", IonoModel::perfectly_modeled(40.0)),
        ("iono_half", IonoModel::klobuchar_like(40.0)),
        (
            "iono_none",
            IonoModel {
                vtec_tecu: 40.0,
                model_error: 1.0,
            },
        ),
    ] {
        let sim = Simulator::new(
            Scenario::ideal_static(name)
                .with_duration(900.0, 30.0)
                .with_iono(iono),
            frame,
        );

        let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Absolute));
        check(
            &report,
            report.stats().rms_3d_m < 0.02,
            "iono free combination",
        );
        assert_contract(Preset::CPP, &report, 0.0);
        rms.push(report.stats().rms_3d_m);
    }

    assert!((rms[0] - rms[2]).abs() < 0.005);
}

/// SPP is sensitive to the ionosphere model, as expected (physics check of the test bench).
#[test]
fn spp_ionosphere_sensitivity() {
    let frame = frame();

    let run_with = |name: &str, iono: IonoModel| {
        let sim = Simulator::new(
            Scenario::ideal_static(name)
                .with_duration(900.0, 30.0)
                .with_iono(iono),
            frame,
        );
        run(
            &sim,
            &RunOptions::new(Method::SPP, Technique::Absolute)
                .with_params(params(UserProfile::Static))
                .with_config(no_postfit),
        )
    };

    let perfect = run_with("spp_iono_perfect", IonoModel::perfectly_modeled(20.0));
    let half = run_with("spp_iono_half", IonoModel::klobuchar_like(20.0));
    let none = run_with(
        "spp_iono_none",
        IonoModel {
            vtec_tecu: 20.0,
            model_error: 1.0,
        },
    );

    check(
        &perfect,
        perfect.stats().rms_3d_m < 0.02,
        "perfect iono model",
    );
    assert_contract(Preset::SPP, &perfect, 0.0);
    // Klobuchar like model (50 % error): nominal SPP conditions
    assert_contract(Preset::SPP, &half, 0.0);

    // unmodeled ionosphere mostly projects onto the vertical (positive) and the clock
    let (h, n) = (half.stats(), none.stats());
    check(
        &half,
        h.mean_enu_m.2 > 1.0 && h.mean_enu_m.2 < 6.0,
        "half iono: vertical bias",
    );
    check(
        &none,
        n.mean_enu_m.2 > 1.5 * h.mean_enu_m.2,
        "no iono: about twice the bias",
    );
    check(&none, n.rms_h_m < 1.0, "horizontal is mostly preserved");
}

/// Unmodeled zenith wet delay biases the vertical (motivates ZWD estimation, PPP milestone).
#[test]
fn unmodeled_troposphere_biases_the_vertical() {
    let frame = frame();

    let sim = Simulator::new(
        Scenario::ideal_static("zwd").with_tropo(TropoModel {
            residual_zwd_m: 0.10,
        }),
        frame,
    );

    let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Absolute));
    let stats = report.stats();

    check(&report, stats.rms_h_m < 0.1, "horizontal");
    check(
        &report,
        stats.mean_enu_m.2 > 0.1 && stats.mean_enu_m.2 < 0.6,
        "10 cm ZWD -> 0.1 .. 0.6 m vertical bias",
    );
}

/// Same scenario, same configuration: bit identical results.
#[test]
fn runs_are_deterministic() {
    let frame = frame();
    let scenario = Scenario::realistic_static("determinism").with_duration(600.0, 30.0);
    let opts = RunOptions::new(Method::CPP, Technique::Absolute);

    let a = run(&Simulator::new(scenario.clone(), frame), &opts);
    let b = run(&Simulator::new(scenario, frame), &opts);

    assert_eq!(a.to_csv(), b.to_csv());
}
