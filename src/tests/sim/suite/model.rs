//! White box checks: simulator self consistency, and solver building blocks
//! (orbit fit, transmission time, relativistic models) against the truth.
use nalgebra::Vector3;

use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    pool::Pool,
    prelude::{Carrier, Config, Constellation, Rc, SV},
    tests::{
        almanac,
        ephemeris::NullEph,
        sim::{
            ClockModel, CycleSlip, IonoModel, NoiseModel, Scenario, SimConstellation, Simulator,
            Station,
        },
        time::NullTime,
    },
};

use super::frame;

fn e(prn: u8) -> SV {
    SV::new(Constellation::Galileo, prn)
}

#[test]
fn simulation_is_reproducible() {
    let frame = frame();
    let scenario = Scenario::realistic_static("repro");

    let a = Simulator::new(scenario.clone(), frame);
    let b = Simulator::new(scenario.clone(), frame);
    let c = Simulator::new(scenario.with_seed(2), frame);

    for k in [0, 10, 100] {
        let (ca, cb, cc) = (
            a.observe(Station::Rover, k),
            b.observe(Station::Rover, k),
            c.observe(Station::Rover, k),
        );

        assert_eq!(ca.len(), cb.len());

        for (x, y) in ca.iter().zip(cb.iter()) {
            assert_eq!(
                format!("{:?}", x.observations),
                format!("{:?}", y.observations)
            );
        }

        // other seed: other noise
        assert_ne!(
            format!("{:?}", ca[0].observations),
            format!("{:?}", cc[0].observations)
        );
    }
}

#[test]
fn constellation_visibility_over_one_day() {
    let frame = frame();

    for (constellation, min_expected) in
        [(SimConstellation::Galileo, 5), (SimConstellation::GPS, 4)]
    {
        let scenario = Scenario::ideal_static("visibility")
            .with_constellations(&[constellation])
            .with_duration(86_400.0, 600.0);

        let sim = Simulator::new(scenario.clone(), frame);

        let mut min_visible = usize::MAX;
        let mut max_visible = 0;

        for k in 0..scenario.num_epochs() {
            // count above the default solver mask (12.5°)
            let visible = sim
                .observe(Station::Rover, k)
                .iter()
                .filter(|cd| {
                    sim.signal_truth(Station::Rover, cd.sv, k)
                        .unwrap()
                        .elevation_deg
                        > 12.5
                })
                .count();

            min_visible = min_visible.min(visible);
            max_visible = max_visible.max(visible);
        }

        assert!(
            min_visible >= min_expected,
            "{constellation:?}: only {min_visible} satellites visible at worst"
        );
        assert!(
            max_visible <= 14,
            "{constellation:?}: {max_visible} visible satellites"
        );
    }
}

/// Noise free observations must match the forward model exactly.
#[test]
fn forward_model_identity() {
    let frame = frame();

    let scenario = Scenario::ideal_static("forward")
        .with_iono(IonoModel::perfectly_modeled(30.0))
        .with_rover_clock(ClockModel::tcxo());

    let sim = Simulator::new(scenario, frame);

    for k in [0, 50, 120] {
        let dt_r = sim.station_clock_s(Station::Rover, k);

        for cd in sim.observe(Station::Rover, k) {
            let truth = sim.signal_truth(Station::Rover, cd.sv, k).unwrap();

            let common = truth.geometric_range_m
                + truth.shapiro_m
                + SPEED_OF_LIGHT_M_S * (dt_r - truth.sat_clock_s)
                + truth.tropo_m;

            for obs in cd.observations.iter() {
                let f_hz = obs.carrier.frequency_hz();
                let lambda = SPEED_OF_LIGHT_M_S / f_hz;
                let iono = sim.iono_delay_m(truth.elevation_deg, f_hz);
                let n = if obs.carrier == Carrier::L1 {
                    truth.ambiguity_1
                } else {
                    truth.ambiguity_2
                };

                let code = obs.pseudo_range_m.unwrap();
                let phase = obs.phase_range_m.unwrap();

                assert!((code - (common + iono)).abs() < 1E-6, "{} code", cd.sv);
                assert!(
                    (phase - (common - iono + lambda * n as f64)).abs() < 1E-6,
                    "{} phase",
                    cd.sv
                );

                // first order ionosphere: f² scaling
                let iono_l1 = sim.iono_delay_m(truth.elevation_deg, Carrier::L1.frequency_hz());
                let expected_ratio = (Carrier::L1.frequency_hz() / f_hz).powi(2);
                assert!((iono / iono_l1 - expected_ratio).abs() < 1E-9);
            }

            // light time consistency
            assert!(
                (truth.time_of_flight_s * SPEED_OF_LIGHT_M_S - truth.geometric_range_m).abs()
                    < 1E-6
            );
        }
    }
}

#[test]
fn noise_statistics() {
    let frame = frame();

    let noise = NoiseModel {
        code_sigma_m: 0.5,
        phase_sigma_m: 0.002,
        elevation_dependent: false,
    };

    let ideal = Simulator::new(
        Scenario::ideal_static("n0").with_duration(7200.0, 30.0),
        frame,
    );
    let noisy = Simulator::new(
        Scenario::ideal_static("n1")
            .with_duration(7200.0, 30.0)
            .with_noise(noise),
        frame,
    );

    let (mut code_sq, mut phase_sq, mut n) = (0.0, 0.0, 0);

    for k in 0..ideal.scenario().num_epochs() {
        let (a, b) = (
            ideal.observe(Station::Rover, k),
            noisy.observe(Station::Rover, k),
        );
        for (x, y) in a.iter().zip(b.iter()) {
            for (ox, oy) in x.observations.iter().zip(y.observations.iter()) {
                code_sq += (oy.pseudo_range_m.unwrap() - ox.pseudo_range_m.unwrap()).powi(2);
                phase_sq += (oy.phase_range_m.unwrap() - ox.phase_range_m.unwrap()).powi(2);
                n += 1;
            }
        }
    }

    let code_rms = (code_sq / n as f64).sqrt();
    let phase_rms = (phase_sq / n as f64).sqrt();

    assert!(n > 2000);
    assert!((code_rms - 0.5).abs() < 0.03, "code rms {code_rms}");
    assert!((phase_rms - 0.002).abs() < 0.0002, "phase rms {phase_rms}");
}

#[test]
fn ambiguities_and_cycle_slips() {
    let frame = frame();

    let scenario = Scenario::ideal_static("slips").with_cycle_slip(CycleSlip {
        sv: e(2),
        at_s: 600.0,
        cycles_1: 7,
        cycles_2: -3,
    });

    let sim = Simulator::new(scenario, frame);

    let before_1 = sim.ambiguity(Station::Rover, e(2), false, 599.0);
    let before_2 = sim.ambiguity(Station::Rover, e(2), true, 599.0);

    assert_eq!(
        sim.ambiguity(Station::Rover, e(2), false, 0.0),
        before_1,
        "ambiguity must be constant"
    );
    assert_eq!(
        sim.ambiguity(Station::Rover, e(2), false, 600.0),
        before_1 + 7
    );
    assert_eq!(
        sim.ambiguity(Station::Rover, e(2), true, 3000.0),
        before_2 - 3
    );

    // other satellites and the base are not affected
    assert_eq!(
        sim.ambiguity(Station::Rover, e(3), false, 0.0),
        sim.ambiguity(Station::Rover, e(3), false, 3000.0)
    );
    assert_ne!(
        sim.ambiguity(Station::Rover, e(2), false, 0.0),
        sim.ambiguity(Station::Base, e(2), false, 0.0)
    );
}

/// The solver orbit fit (transmission time + Sagnac rotation) must match the
/// simulated truth, whatever the receiver clock offset (once it is known).
#[test]
fn pool_orbit_fit_matches_truth() {
    let frame = frame();
    let dt_r = 5.0E-4;

    let sim = Simulator::new(
        Scenario::ideal_static("orbits").with_rover_clock(ClockModel {
            offset_s: dt_r,
            drift_s_s: 0.0,
            random_walk_s_sqrt_s: 0.0,
        }),
        frame,
    );

    let mut pool = Pool::allocate(
        almanac(),
        Config::default(),
        frame,
        Rc::new(NullEph {}),
        Rc::new(sim.orbit_source()),
        Rc::new(sim.environment()),
        Rc::new(sim.spaceborn_biases()),
    );

    for (clock_estimate_s, tolerance_m, must_fail) in [(dt_r, 1.0E-3, false), (0.0, 1.0E-3, true)] {
        let k = 10;
        pool.new_epoch(&sim.observe(Station::Rover, k));
        pool.pre_fit("rover", &NullTime {});
        pool.orbital_states_fit_with_clock("rover", clock_estimate_s);

        assert!(pool.len() > 4);

        let mut max_err: f64 = 0.0;

        for cd in pool.candidates() {
            let truth = sim.signal_truth(Station::Rover, cd.sv, k).unwrap();
            let fitted = cd.orbit.unwrap().to_cartesian_pos_vel() * 1.0E3;
            let err = (Vector3::new(fitted[0], fitted[1], fitted[2]) - truth.sat_position_m).norm();
            max_err = max_err.max(err);
        }

        if must_fail {
            // ignoring a 0.5 ms receiver clock offset misplaces satellites by ~1 m
            assert!(max_err > 0.3, "the receiver clock should matter: {max_err}");
        } else {
            assert!(max_err < tolerance_m, "orbit fit error {max_err} m");
        }
    }
}

/// Solver Shapiro model against the (independent) simulator formula
#[test]
fn shapiro_model_matches_truth() {
    use crate::pool::postfit::shapiro_delay_m;

    let frame = frame();
    let sim = Simulator::new(Scenario::ideal_static("shapiro"), frame);
    let rx = sim.station_position_m(Station::Rover, 0.0);

    for cd in sim.observe(Station::Rover, 0) {
        let truth = sim.signal_truth(Station::Rover, cd.sv, 0).unwrap();
        let model = shapiro_delay_m(
            truth.sat_position_m.norm(),
            rx.norm(),
            truth.geometric_range_m,
        );
        assert!((model - truth.shapiro_m).abs() < 1E-9);
        assert!(truth.shapiro_m > 0.01 && truth.shapiro_m < 0.025);
    }
}
