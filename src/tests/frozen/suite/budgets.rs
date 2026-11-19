//! Performance budgets on the frozen real data.
//!
//! Every (dataset, technique, constellations) case is replayed and checked against
//! its budget: availability, horizontal / vertical RMS and maximal 3D error,
//! statistics computed after the initial convergence (10' dense, 1 h day).
//! Budgets are max(1.4 x, +5 cm) the values observed when they were set: the
//! processing is deterministic, any degradation fails. Tighten them when the
//! solver improves.
//!
//! Every case must also comply with the contract of its preset (`tests::spec`,
//! published in the README), unless a derogation (known limitation) is listed.
//!
//! Known limitations are expressed as ignored `target_*` tests.
use crate::{
    prelude::{
        Constellation::{self, BeiDou, GPS, Galileo},
        Method,
    },
    tests::{
        frozen::{Dataset, Replay, replay},
        spec::Preset,
    },
};

struct Budget {
    dataset: Dataset,
    method: Method,
    rtk: bool,
    constellations: &'static [Constellation],
    min_availability: f64,
    max_rms_h_m: f64,
    max_rms_v_m: f64,
    max_3d_m: f64,
    /// Known violation of the preset contract (reason): the contract is then
    /// checked by the target_frozen_contract_derogations test only.
    derogation: Option<&'static str>,
}

#[rustfmt::skip]
const BUDGETS: &[Budget] = &[
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 2.39, max_rms_v_m: 7.68, max_3d_m: 8.93, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[Galileo], min_availability: 1.00, max_rms_h_m: 16.66, max_rms_v_m: 9.88, max_3d_m: 22.62, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[BeiDou], min_availability: 0.73, max_rms_h_m: 1.09, max_rms_v_m: 2.05, max_3d_m: 2.43, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 7.10, max_rms_v_m: 1.89, max_3d_m: 9.29, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 4.84, max_rms_v_m: 1.85, max_3d_m: 6.80, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 1.30, max_rms_v_m: 1.14, max_3d_m: 1.95, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[Galileo], min_availability: 1.00, max_rms_h_m: 0.80, max_rms_v_m: 0.19, max_3d_m: 0.97, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 1.26, max_rms_v_m: 1.21, max_3d_m: 2.51, derogation: Some("GPS / Galileo inter-system bias not estimated (single receiver clock)") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.29, max_rms_v_m: 0.14, max_3d_m: 0.36, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[Galileo], min_availability: 0.87, max_rms_h_m: 0.45, max_rms_v_m: 0.28, max_3d_m: 0.60, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.33, max_rms_v_m: 0.31, max_3d_m: 0.54, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.35, max_rms_v_m: 0.29, max_3d_m: 0.67, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[Galileo], min_availability: 0.87, max_rms_h_m: 0.65, max_rms_v_m: 0.10, max_3d_m: 0.70, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.49, max_rms_v_m: 0.25, max_3d_m: 0.69, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 2.29, max_rms_v_m: 2.51, max_3d_m: 6.45, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[Galileo], min_availability: 0.91, max_rms_h_m: 4.14, max_rms_v_m: 4.01, max_3d_m: 13.20, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[BeiDou], min_availability: 0.43, max_rms_h_m: 1.26, max_rms_v_m: 3.17, max_3d_m: 3.98, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 2.08, max_rms_v_m: 1.83, max_3d_m: 6.17, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 1.44, max_rms_v_m: 2.45, max_3d_m: 4.32, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 1.22, max_rms_v_m: 0.45, max_3d_m: 1.90, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[Galileo], min_availability: 0.91, max_rms_h_m: 1.08, max_rms_v_m: 0.24, max_3d_m: 1.17, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 1.36, max_rms_v_m: 0.74, max_3d_m: 2.56, derogation: Some("GPS / Galileo inter-system bias not estimated (single receiver clock)") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[GPS], min_availability: 0.95, max_rms_h_m: 0.33, max_rms_v_m: 0.11, max_3d_m: 0.43, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[Galileo], min_availability: 0.64, max_rms_h_m: 0.32, max_rms_v_m: 0.17, max_3d_m: 0.43, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[BeiDou], min_availability: 0.21, max_rms_h_m: 0.47, max_rms_v_m: 0.78, max_3d_m: 2.02, derogation: Some("BeiDou: broadcast orbits and clocks, 4 to 6 satellites") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.32, max_rms_v_m: 0.14, max_3d_m: 0.50, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[GPS], min_availability: 0.95, max_rms_h_m: 0.37, max_rms_v_m: 0.17, max_3d_m: 0.70, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[Galileo], min_availability: 0.64, max_rms_h_m: 0.38, max_rms_v_m: 0.15, max_3d_m: 0.79, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[BeiDou], min_availability: 0.21, max_rms_h_m: 0.79, max_rms_v_m: 1.88, max_3d_m: 5.40, derogation: Some("BeiDou: broadcast orbits and clocks, 4 to 6 satellites") },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.41, max_rms_v_m: 0.13, max_3d_m: 0.86, derogation: None },
];

fn convergence_s(dataset: Dataset) -> f64 {
    match dataset {
        Dataset::Dense => 600.0,
        Dataset::Day => 3600.0,
    }
}

fn check_budgets(dataset: Dataset, method: Method, rtk: bool) {
    crate::tests::init_logger();

    let mut failures = Vec::new();
    let mut count = 0;

    for budget in BUDGETS
        .iter()
        .filter(|b| b.dataset == dataset && b.method == method && b.rtk == rtk)
    {
        count += 1;

        let report = replay(&Replay::new(dataset, method, rtk, budget.constellations));
        let stats = report.stats_from(convergence_s(dataset));

        let checks = [
            (
                "availability",
                report.availability(),
                budget.min_availability,
                report.availability() >= budget.min_availability,
            ),
            (
                "rms horizontal (m)",
                stats.rms_h_m,
                budget.max_rms_h_m,
                stats.rms_h_m <= budget.max_rms_h_m,
            ),
            (
                "rms vertical (m)",
                stats.rms_u_m,
                budget.max_rms_v_m,
                stats.rms_u_m <= budget.max_rms_v_m,
            ),
            (
                "max 3D (m)",
                stats.max_3d_m,
                budget.max_3d_m,
                stats.max_3d_m <= budget.max_3d_m,
            ),
        ];

        // preset contract (README), unless a known derogation applies
        if budget.derogation.is_none() {
            let preset = Preset::new(method, rtk);
            for violation in preset
                .contract()
                .report_violations(&report, convergence_s(dataset))
                .into_iter()
                .take(3)
            {
                failures.push(format!("{}: {preset:?} contract: {violation}", report.name));
            }
        }

        for (what, value, limit, ok) in checks {
            if !ok || !value.is_finite() {
                failures.push(format!(
                    "{}: {what} = {value:.3} (budget {limit:.3})",
                    report.name
                ));
            }
        }
    }

    assert!(count > 0, "no budget defined");
    assert!(
        failures.is_empty(),
        "budgets exceeded:\n{}",
        failures.join("\n")
    );
}

#[test]
fn frozen_dense_spp() {
    check_budgets(Dataset::Dense, Method::SPP, false);
}

#[test]
fn frozen_dense_cpp() {
    check_budgets(Dataset::Dense, Method::CPP, false);
}

#[test]
fn frozen_dense_rtk_spp() {
    check_budgets(Dataset::Dense, Method::SPP, true);
}

#[test]
fn frozen_dense_rtk_cpp() {
    check_budgets(Dataset::Dense, Method::CPP, true);
}

#[test]
fn frozen_day_spp() {
    check_budgets(Dataset::Day, Method::SPP, false);
}

#[test]
fn frozen_day_cpp() {
    check_budgets(Dataset::Day, Method::CPP, false);
}

#[test]
fn frozen_day_rtk_spp() {
    check_budgets(Dataset::Day, Method::SPP, true);
}

#[test]
fn frozen_day_rtk_cpp() {
    check_budgets(Dataset::Day, Method::CPP, true);
}

/// Frozen data sanity: every epoch has data for the three constellations,
/// orbits and clocks exist for every observed satellite.
#[test]
fn frozen_data_consistency() {
    for dataset in [Dataset::Dense, Dataset::Day] {
        let data = dataset.load();

        assert!(
            data.rover.epochs.len() > 100,
            "{dataset:?}: not enough epochs"
        );
        assert_eq!(data.rover.name, "ESBC00DNK");
        assert_eq!(data.base.name, "MOJN00DNK");

        // 65 km baseline
        let baseline_m = (data.rover.arp_ecef_m() - data.base.arp_ecef_m()).norm();
        assert!(
            (baseline_m - 65_000.0).abs() < 1_000.0,
            "baseline {baseline_m} m"
        );

        for (epoch, candidates) in data.rover.epochs.iter() {
            for constellation in [GPS, Galileo, BeiDou] {
                assert!(
                    candidates
                        .iter()
                        .any(|cd| cd.sv.constellation == constellation),
                    "{dataset:?} {epoch}: no {constellation} data"
                );
            }

            for cd in candidates {
                assert!(
                    data.clocks.get(cd.sv, *epoch).is_some(),
                    "{dataset:?} {epoch}: missing clock for {}",
                    cd.sv
                );
                assert!(
                    data.orbits.position_km(cd.sv, *epoch).is_some(),
                    "{dataset:?} {epoch}: missing orbit for {}",
                    cd.sv
                );
            }
        }
    }
}

/// Target: lift the contract derogations. Each derogated case must eventually
/// comply with its preset contract (see the derogation reasons in [BUDGETS]).
#[test]
#[ignore]
fn target_frozen_contract_derogations() {
    let mut failures = Vec::new();

    for budget in BUDGETS.iter().filter(|b| b.derogation.is_some()) {
        let report = replay(&Replay::new(
            budget.dataset,
            budget.method,
            budget.rtk,
            budget.constellations,
        ));

        let violations = Preset::new(budget.method, budget.rtk)
            .contract()
            .report_violations(&report, convergence_s(budget.dataset));

        if !violations.is_empty() {
            failures.push(format!(
                "{} ({}): {} violations, first: {}",
                report.name,
                budget.derogation.unwrap(),
                violations.len(),
                violations[0]
            ));
        }
    }

    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Target: BeiDou CPP. The B1I group delay (TGD1) is applied as is to the iono free
/// combination (broadcast clocks refer to B3I): ~25 m vertical bias today.
#[test]
#[ignore]
fn target_frozen_cpp_beidou() {
    let report = replay(&Replay::new(Dataset::Day, Method::CPP, false, &[BeiDou]));
    let stats = report.stats_from(3600.0);
    assert!(stats.rms_3d_m < 2.0, "{report}");
}

/// Target: multi constellation RTK. A single pivot satellite for all constellations
/// (different signals) leads to failures (availability 86 %).
#[test]
#[ignore]
fn target_frozen_rtk_multi_constellation() {
    for method in [Method::SPP, Method::CPP] {
        let report = replay(&Replay::new(
            Dataset::Day,
            method,
            true,
            &[GPS, Galileo, BeiDou],
        ));
        assert!(report.availability() > 0.99, "{report}");
    }
}

/// Target: SPP with an ionosphere model (Klobuchar from the broadcast message).
/// The test environment does not model the ionosphere today.
#[test]
#[ignore]
fn target_frozen_spp_with_ionosphere_model() {
    let report = replay(&Replay::new(Dataset::Day, Method::SPP, false, &[GPS]));
    assert!(report.stats_from(3600.0).rms_3d_m < 1.5, "{report}");
}

/// Prints a markdown table of every frozen case (CI job summary).
#[test]
#[ignore]
fn frozen_performance_report() {
    let mut table = String::from(
        "| case | avail | rms H (m) | rms V (m) | rms 3D (m) | p95 3D (m) | max 3D (m) | budget H / V / max (m) | contract |\n|---|---|---|---|---|---|---|---|---|\n",
    );

    for budget in BUDGETS {
        let report = replay(&Replay::new(
            budget.dataset,
            budget.method,
            budget.rtk,
            budget.constellations,
        ));
        let s = report.stats_from(convergence_s(budget.dataset));

        let violations = Preset::new(budget.method, budget.rtk)
            .contract()
            .report_violations(&report, convergence_s(budget.dataset));

        let contract = match (violations.is_empty(), budget.derogation) {
            (true, _) => "OK".to_string(),
            (false, Some(reason)) => format!("derogation: {reason}"),
            (false, None) => format!("**{} violations**", violations.len()),
        };

        table.push_str(&format!(
            "| {} | {:.1}% | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.2} / {:.2} / {:.2} | {} |\n",
            report.name,
            report.availability() * 100.0,
            s.rms_h_m,
            s.rms_u_m,
            s.rms_3d_m,
            s.p95_3d_m,
            s.max_3d_m,
            budget.max_rms_h_m,
            budget.max_rms_v_m,
            budget.max_3d_m,
            contract,
        ));
    }

    println!("\n{table}");

    if let Ok(dir) = std::env::var("GNSS_RTK_SIM_OUTPUT") {
        let _ = std::fs::create_dir_all(&dir);
        let _ = std::fs::write(std::path::Path::new(&dir).join("frozen_report.md"), &table);
    }
}
