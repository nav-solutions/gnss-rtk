//! Performance budgets on the frozen real data.
//!
//! Every (dataset, technique, constellations) case is replayed and checked against
//! its budget: availability, horizontal / vertical RMS and maximal 3D error,
//! statistics computed after the initial convergence (10' dense, 1 h day).
//! Truth: ITRF2014 for the absolute techniques, ETRF2000 for RTK
//! (see [crate::tests::frozen::StationData::itrf2014_arp_ecef_m]).
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
        spec::{PPP_STEADY_STATE_S, Preset},
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
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 3.18, max_rms_v_m: 7.37, max_3d_m: 8.72, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[Galileo], min_availability: 1.00, max_rms_h_m: 16.24, max_rms_v_m: 10.39, max_3d_m: 22.73, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[BeiDou], min_availability: 0.75, max_rms_h_m: 0.35, max_rms_v_m: 2.61, max_3d_m: 2.75, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 6.81, max_rms_v_m: 2.06, max_3d_m: 9.20, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 4.51, max_rms_v_m: 2.15, max_3d_m: 6.72, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.38, max_rms_v_m: 1.43, max_3d_m: 1.58, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[Galileo], min_availability: 1.00, max_rms_h_m: 0.14, max_rms_v_m: 0.48, max_3d_m: 0.78, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 1.29, max_rms_v_m: 1.57, max_3d_m: 3.10, derogation: Some("GPS / Galileo inter-system bias not estimated (single receiver clock)") },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.29, max_rms_v_m: 0.14, max_3d_m: 0.36, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[Galileo], min_availability: 0.87, max_rms_h_m: 0.45, max_rms_v_m: 0.28, max_3d_m: 0.60, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::SPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.33, max_rms_v_m: 0.31, max_3d_m: 0.54, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.35, max_rms_v_m: 0.29, max_3d_m: 0.67, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[Galileo], min_availability: 0.87, max_rms_h_m: 0.65, max_rms_v_m: 0.10, max_3d_m: 0.70, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::CPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.49, max_rms_v_m: 0.25, max_3d_m: 0.69, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 1.98, max_rms_v_m: 2.44, max_3d_m: 6.83, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[Galileo], min_availability: 0.93, max_rms_h_m: 4.52, max_rms_v_m: 3.97, max_3d_m: 13.42, derogation: Some("Galileo E24: ~12 m E1 code bias not rejected (no RAIM), ionosphere not modeled") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[BeiDou], min_availability: 0.45, max_rms_h_m: 0.26, max_rms_v_m: 3.03, max_3d_m: 3.66, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 2.01, max_rms_v_m: 1.78, max_3d_m: 6.05, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 1.56, max_rms_v_m: 2.35, max_3d_m: 4.24, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.27, max_rms_v_m: 0.39, max_3d_m: 1.42, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[Galileo], min_availability: 0.93, max_rms_h_m: 0.11, max_rms_v_m: 0.17, max_3d_m: 0.32, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.57, max_rms_v_m: 0.73, max_3d_m: 1.90, derogation: Some("GPS / Galileo inter-system bias not estimated (single receiver clock)") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[GPS], min_availability: 0.95, max_rms_h_m: 0.33, max_rms_v_m: 0.11, max_3d_m: 0.43, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[Galileo], min_availability: 0.64, max_rms_h_m: 0.32, max_rms_v_m: 0.17, max_3d_m: 0.43, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[BeiDou], min_availability: 0.21, max_rms_h_m: 0.47, max_rms_v_m: 0.78, max_3d_m: 2.02, derogation: Some("BeiDou: broadcast orbits and clocks, 4 to 6 satellites") },
    Budget { dataset: Dataset::Day, method: Method::SPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.32, max_rms_v_m: 0.14, max_3d_m: 0.50, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[GPS], min_availability: 0.95, max_rms_h_m: 0.37, max_rms_v_m: 0.17, max_3d_m: 0.70, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[Galileo], min_availability: 0.64, max_rms_h_m: 0.38, max_rms_v_m: 0.15, max_3d_m: 0.79, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[BeiDou], min_availability: 0.21, max_rms_h_m: 0.79, max_rms_v_m: 1.88, max_3d_m: 5.40, derogation: Some("BeiDou: broadcast orbits and clocks, 4 to 6 satellites") },
    Budget { dataset: Dataset::Day, method: Method::CPP, rtk: true, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.41, max_rms_v_m: 0.13, max_3d_m: 0.86, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::PPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.20, max_rms_v_m: 0.14, max_3d_m: 0.32, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::PPP, rtk: false, constellations: &[Galileo], min_availability: 1.00, max_rms_h_m: 0.15, max_rms_v_m: 0.11, max_3d_m: 0.23, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::PPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.13, max_rms_v_m: 0.10, max_3d_m: 0.17, derogation: None },
    Budget { dataset: Dataset::Dense, method: Method::PPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 0.34, max_rms_v_m: 0.22, max_3d_m: 0.65, derogation: Some("BeiDou: broadcast orbits and clocks (decimeter level) mixed with precise products") },
    Budget { dataset: Dataset::Day, method: Method::PPP, rtk: false, constellations: &[GPS], min_availability: 1.00, max_rms_h_m: 0.09, max_rms_v_m: 0.14, max_3d_m: 0.22, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::PPP, rtk: false, constellations: &[Galileo], min_availability: 0.90, max_rms_h_m: 0.10, max_rms_v_m: 0.13, max_3d_m: 0.20, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::PPP, rtk: false, constellations: &[GPS, Galileo], min_availability: 1.00, max_rms_h_m: 0.10, max_rms_v_m: 0.14, max_3d_m: 0.21, derogation: None },
    Budget { dataset: Dataset::Day, method: Method::PPP, rtk: false, constellations: &[GPS, Galileo, BeiDou], min_availability: 1.00, max_rms_h_m: 0.12, max_rms_v_m: 0.18, max_3d_m: 0.41, derogation: None },
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

            // steady state (converging presets), on the 24 h dataset. Scope: precise
            // products only (BeiDou: broadcast orbits and clocks in these data)
            let precise_products = !budget.constellations.contains(&BeiDou);

            if let (Dataset::Day, Some(steady), true) =
                (dataset, preset.steady_state_contract(), precise_products)
            {
                for violation in steady
                    .report_violations(&report, PPP_STEADY_STATE_S)
                    .into_iter()
                    .take(3)
                {
                    failures.push(format!(
                        "{}: {preset:?} steady state contract: {violation}",
                        report.name
                    ));
                }
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

#[test]
fn frozen_dense_ppp() {
    check_budgets(Dataset::Dense, Method::PPP, false);
}

#[test]
fn frozen_day_ppp() {
    check_budgets(Dataset::Day, Method::PPP, false);
}

/// PPP receiver clock precision on real data, without clock truth: two
/// independent estimates (GPS only, Galileo only: different satellites, orbits
/// and clocks) must agree, up to a constant inter-system bias, at the phase level.
/// Their common behavior (~1 ns white noise, ~5 ns wander per day) is the receiver
/// oscillator itself.
#[test]
fn frozen_ppp_clock_consistency() {
    crate::tests::init_logger();

    // (dataset, max std, max deviation) in ns: over a day, the inter-system bias
    // itself varies (clock products datum, receiver temperature)
    for (dataset, max_std, max_dev) in [(Dataset::Dense, 0.2, 0.6), (Dataset::Day, 0.6, 1.2)] {
        let gps = replay(&Replay::new(dataset, Method::PPP, false, &[GPS]));
        let gal = replay(&Replay::new(dataset, Method::PPP, false, &[Galileo]));

        let differences_ns = gps
            .outcomes
            .iter()
            .zip(gal.outcomes.iter())
            .filter(|(a, _)| a.t_s >= convergence_s(dataset))
            .filter_map(|(a, b)| match (a.solution, b.solution) {
                (Some((_, dt_a)), Some((_, dt_b))) => Some((dt_a - dt_b) * 1.0E9),
                _ => None,
            })
            .collect::<Vec<_>>();

        assert!(
            differences_ns.len() > 100,
            "{dataset:?}: not enough solutions"
        );

        let n = differences_ns.len() as f64;
        let mean = differences_ns.iter().sum::<f64>() / n;
        let std = (differences_ns
            .iter()
            .map(|d| (d - mean).powi(2))
            .sum::<f64>()
            / n)
            .sqrt();
        let max = differences_ns
            .iter()
            .map(|d| (d - mean).abs())
            .fold(0.0, f64::max);

        println!(
            "{dataset:?}: GPS - Galileo clock: bias {mean:.3} ns, std {std:.3} ns, max {max:.3} ns"
        );

        assert!(
            std < max_std,
            "{dataset:?}: GPS / Galileo PPP clocks std {std:.3} ns"
        );
        assert!(
            max < max_dev,
            "{dataset:?}: GPS / Galileo PPP clocks max deviation {max:.3} ns"
        );
    }
}

/// The RINEX header coordinates (ETRF2000) and the ITRF2014 truth of the
/// absolute techniques: plate motion since 1989 (~2.5 cm/yr towards the
/// North-East in Denmark), the baseline (RTK) is preserved.
#[test]
fn frozen_reference_frames() {
    use crate::tests::sim::geodesy::ecef_delta_to_enu;

    for dataset in [Dataset::Dense, Dataset::Day] {
        let data = dataset.load();

        for station in [&data.rover, &data.base] {
            let etrf = station.arp_ecef_m();
            let itrf = station.itrf2014_arp_ecef_m();
            let enu = ecef_delta_to_enu(&etrf, &(itrf - etrf));

            assert!(
                (enu.x - 0.52).abs() < 0.02,
                "{}: east shift {}",
                station.name,
                enu.x
            );
            assert!(
                (enu.y - 0.56).abs() < 0.02,
                "{}: north shift {}",
                station.name,
                enu.y
            );
            assert!(enu.z.abs() < 0.03, "{}: up shift {}", station.name, enu.z);
        }

        let baseline_etrf = data.rover.arp_ecef_m() - data.base.arp_ecef_m();
        let baseline_itrf = data.rover.itrf2014_arp_ecef_m() - data.base.itrf2014_arp_ecef_m();
        // frame rotation (~20 mas) over 65 km: a few mm
        assert!((baseline_etrf - baseline_itrf).norm() < 0.01);
    }
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
