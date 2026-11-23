//! Performance report: runs a catalogue of scenarios and prints a markdown table.
//!
//! ```bash
//! cargo test --release performance_report -- --ignored --nocapture
//! # optional: CSV per run + report.md
//! GNSS_RTK_SIM_OUTPUT=target/sim cargo test --release performance_report -- --ignored
//! ```
//!
//! Track this report across releases to monitor the performance of each technique.
use crate::{
    prelude::{Method, UserProfile},
    tests::sim::{
        IonoModel, RunOptions, Scenario, SimConstellation, Simulator, Technique, Trajectory,
        geodesy::geodetic_to_ecef, run, scenario::REFERENCE_SITE_LAT_LON_H,
    },
};

use super::{frame, no_postfit, params};

struct Entry {
    scenario: Scenario,
    opts: RunOptions,
    /// statistics are computed after this instant (s)
    from_s: f64,
}

fn catalogue() -> Vec<Entry> {
    let (lat, lon, h) = REFERENCE_SITE_LAT_LON_H;
    let center = geodetic_to_ecef(lat.to_radians(), lon.to_radians(), h);

    let static_opts = |method: Method, technique: Technique| {
        RunOptions::new(method, technique)
            .with_params(params(UserProfile::Static))
            .with_config(no_postfit)
    };

    let mut entries = Vec::new();

    for (method, technique) in [
        (Method::SPP, Technique::Absolute),
        (Method::CPP, Technique::Absolute),
        (Method::SPP, Technique::Rtk),
        (Method::CPP, Technique::Rtk),
    ] {
        let tag = format!(
            "{}-{}",
            if technique == Technique::Rtk {
                "rtk"
            } else {
                "abs"
            },
            method
        )
        .to_lowercase();

        entries.push(Entry {
            scenario: Scenario::ideal_static(&format!("{tag}_ideal_5km"))
                .with_base_offset_enu(3_000.0, 4_000.0, 0.0),
            opts: RunOptions::new(method, technique),
            from_s: 0.0,
        });

        entries.push(Entry {
            scenario: Scenario::realistic_static(&format!("{tag}_static_geodetic_5km"))
                .with_base_offset_enu(3_000.0, 4_000.0, 0.0),
            opts: static_opts(method, technique),
            from_s: 600.0,
        });

        entries.push(Entry {
            scenario: Scenario::realistic_static(&format!("{tag}_static_gps_gal_5km"))
                .with_constellations(&[SimConstellation::GPS, SimConstellation::Galileo])
                .with_base_offset_enu(3_000.0, 4_000.0, 0.0),
            opts: static_opts(method, technique),
            from_s: 600.0,
        });

        entries.push(Entry {
            scenario: Scenario::realistic_static(&format!("{tag}_pedestrian_geodetic_5km"))
                .with_duration(600.0, 1.0)
                .with_trajectory(Trajectory::Circle {
                    center_ecef_m: center,
                    radius_m: 100.0,
                    period_s: 300.0,
                })
                .with_base_offset_enu(3_000.0, 4_000.0, 0.0),
            opts: RunOptions::new(method, technique)
                .with_params(params(UserProfile::Airplane))
                .with_config(no_postfit),
            from_s: 0.0,
        });
    }

    entries.push(Entry {
        scenario: Scenario::realistic_static("abs-spp_static_no_iono_model").with_iono(IonoModel {
            vtec_tecu: 20.0,
            model_error: 1.0,
        }),
        opts: static_opts(Method::SPP, Technique::Absolute),
        from_s: 600.0,
    });

    entries.push(Entry {
        scenario: Scenario::ideal_static("rtk-cpp_ideal_30km")
            .with_base_offset_enu(18_000.0, 24_000.0, 0.0),
        opts: RunOptions::new(Method::CPP, Technique::Rtk),
        from_s: 0.0,
    });

    entries
}

#[test]
#[ignore]
fn performance_report() {
    let frame = frame();

    let mut table = String::from(
        "| scenario | avail | rms E | rms N | rms U | rms 3D | p95 3D | max 3D | conv<1m | clk rms |\n|---|---|---|---|---|---|---|---|---|---|\n",
    );

    for entry in catalogue() {
        let sim = Simulator::new(entry.scenario, frame);
        let report = run(&sim, &entry.opts);
        let s = report.stats_from(entry.from_s);

        let clock = if entry.opts.technique == Technique::Rtk {
            "-".to_string()
        } else {
            format!("{:.2} ns", report.clock_rms_ns_from(entry.from_s))
        };

        table.push_str(&format!(
            "| {} | {:.1}% | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {} | {} |\n",
            report.name,
            report.availability() * 100.0,
            s.rms_e_m,
            s.rms_n_m,
            s.rms_u_m,
            s.rms_3d_m,
            s.p95_3d_m,
            s.max_3d_m,
            report
                .convergence_time_s(1.0)
                .map(|t| format!("{t:.0} s"))
                .unwrap_or("-".to_string()),
            clock,
        ));
    }

    println!("\n{table}");

    if let Ok(dir) = std::env::var("GNSS_RTK_SIM_OUTPUT") {
        let _ = std::fs::create_dir_all(&dir);
        let _ = std::fs::write(std::path::Path::new(&dir).join("report.md"), &table);
    }
}
