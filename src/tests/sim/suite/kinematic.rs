//! Moving receivers.
//!
//! The navigation filter state is position (+clock) only: dynamics are absorbed
//! by the process noise. A high dynamics profile therefore tracks any trajectory
//! (epoch wise solution), while low dynamics profiles lag behind. Tests marked
//! `#[ignore]` are the targets of a position-velocity state (next milestone).
use crate::{
    prelude::{Method, UserProfile},
    tests::sim::{
        NoiseModel, RunOptions, Scenario, Simulator, Technique, Trajectory,
        geodesy::geodetic_to_ecef, run, scenario::REFERENCE_SITE_LAT_LON_H,
    },
};

use nalgebra::Vector3;

use super::{check, frame, no_postfit, params};

fn center() -> Vector3<f64> {
    let (lat, lon, h) = REFERENCE_SITE_LAT_LON_H;
    geodetic_to_ecef(lat.to_radians(), lon.to_radians(), h)
}

/// Pedestrian: 100 m radius circle in 5 minutes (~2 m/s), 1 Hz
fn pedestrian(name: &str) -> Scenario {
    Scenario::ideal_static(name)
        .with_duration(300.0, 1.0)
        .with_trajectory(Trajectory::Circle {
            center_ecef_m: center(),
            radius_m: 100.0,
            period_s: 300.0,
        })
}

/// Car: 500 m radius circle in 150 s (~21 m/s), 1 Hz
fn car(name: &str) -> Scenario {
    Scenario::ideal_static(name)
        .with_duration(300.0, 1.0)
        .with_trajectory(Trajectory::Circle {
            center_ecef_m: center(),
            radius_m: 500.0,
            period_s: 150.0,
        })
}

#[test]
fn trajectory_generator() {
    let traj = Trajectory::Circle {
        center_ecef_m: center(),
        radius_m: 100.0,
        period_s: 300.0,
    };

    for t in [0.0, 37.0, 150.0, 299.0] {
        let d = (traj.position_m(t) - center()).norm();
        assert!((d - 100.0).abs() < 1E-6);
    }

    let speed = (traj.position_m(1.0) - traj.position_m(0.0)).norm();
    assert!((speed - 2.0 * std::f64::consts::PI * 100.0 / 300.0).abs() < 0.01);

    let linear = Trajectory::Linear {
        start_ecef_m: center(),
        velocity_enu_m_s: Vector3::new(3.0, 4.0, 0.0),
    };
    assert!(((linear.position_m(10.0) - center()).norm() - 50.0).abs() < 1E-6);
}

#[test]
fn high_dynamics_profile_tracks_pedestrian_and_car() {
    let frame = frame();

    for (scenario, max_rms_m) in [
        (pedestrian("kin_pedestrian_airplane"), 0.02),
        (car("kin_car_airplane"), 0.15),
    ] {
        let sim = Simulator::new(scenario, frame);

        let report = run(
            &sim,
            &RunOptions::new(Method::CPP, Technique::Absolute)
                .with_params(params(UserProfile::Airplane))
                .with_config(no_postfit),
        );

        check(&report, report.availability() == 1.0, "availability");
        check(&report, report.stats().rms_3d_m < max_rms_m, "3D RMS");
    }
}

#[test]
fn noisy_pedestrian_with_high_dynamics_profile() {
    let frame = frame();

    let sim = Simulator::new(
        pedestrian("kin_pedestrian_noisy").with_noise(NoiseModel::geodetic()),
        frame,
    );

    let report = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Airplane))
            .with_config(no_postfit),
    );

    // epoch wise solution: code noise only
    check(&report, report.stats().rms_3d_m < 6.0, "3D RMS");
    check(&report, report.stats().rms_h_m < 3.0, "horizontal RMS");
}

/// Target (position-velocity state): the pedestrian profile must track a pedestrian.
#[test]
#[ignore]
fn target_pedestrian_profile_tracks_pedestrian() {
    let frame = frame();
    let sim = Simulator::new(pedestrian("target_kin_pedestrian"), frame);

    let report = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Pedestrian))
            .with_config(no_postfit),
    );

    check(
        &report,
        report.stats().rms_3d_m < 0.1,
        "pedestrian profile must track",
    );
}

/// Target (position-velocity state): the car profile must track a car.
#[test]
#[ignore]
fn target_car_profile_tracks_car() {
    let frame = frame();
    let sim = Simulator::new(car("target_kin_car"), frame);

    let report = run(
        &sim,
        &RunOptions::new(Method::CPP, Technique::Absolute)
            .with_params(params(UserProfile::Car))
            .with_config(no_postfit),
    );

    check(
        &report,
        report.stats().rms_3d_m < 0.2,
        "car profile must track",
    );
}

/// Target: the default configuration (postfit denoising enabled) should not
/// prevent a pedestrian from being tracked. Today the postfit smoother
/// (q = 1E-6 m²) lags by tens of meters.
#[test]
#[ignore]
fn target_default_config_tracks_pedestrian() {
    let frame = frame();
    let sim = Simulator::new(pedestrian("target_kin_default"), frame);

    let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Absolute));

    check(
        &report,
        report.stats().rms_3d_m < 1.0,
        "default config must track a pedestrian",
    );
}
