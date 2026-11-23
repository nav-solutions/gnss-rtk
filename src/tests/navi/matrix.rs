use std::str::FromStr;

use anise::{constants::frames::EARTH_J2000, math::Vector6};
use hifitime::Duration;

use nalgebra::U4;

use crate::{
    cfg::Modeling,
    constants::{EARTH_GRAVITATION, SPEED_OF_LIGHT_M_S},
    navigation::{apriori::Apriori, Navigation, state::State},
    prelude::{
        Almanac, Candidate, Carrier, ClockCorrection, Config, Epoch, Error, Method, Observation,
        Orbit, Vector3,
    },
    tests::{
        bias::NullBias,
        gps::{G01, G02, G03, G04},
        time::NullTime,
        REFERENCE_COORDS_ECEF_M,
    },
    time::AbsoluteTime,
};

#[test]
#[cfg(feature = "embed_ephem")]
fn pvt_failures() {
    let cfg = Config::default().with_modeling(Modeling::no_modeling());

    let almanac = crate::tests::almanac();
    let frame = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let coords_ecef_m = Vector3::new(
        REFERENCE_COORDS_ECEF_M.0,
        REFERENCE_COORDS_ECEF_M.1,
        REFERENCE_COORDS_ECEF_M.2,
    );

    let apriori = Apriori::from_ecef_m(coords_ecef_m, Default::default(), frame).unwrap();

    let t = Epoch::from_str("2020-01-01T00:00:00 GPST").unwrap();

    let candidates = vec![Candidate::new(
        G01,
        t,
        vec![Observation::pseudo_range(Carrier::L1, 1.0, None)],
    )];

    let null_bias = NullBias {};

    let null_time = NullTime {};
    let absolute_time = AbsoluteTime::new(null_time);

    let candidates = vec![
        Candidate::new(
            G01,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 1.0, None)],
        ),
        Candidate::new(
            G02,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 2.0, None)],
        ),
    ];

    let candidates = vec![
        Candidate::new(
            G01,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 1.0, None)],
        ),
        Candidate::new(
            G02,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 2.0, None)],
        ),
        Candidate::new(
            G03,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 3.0, None)],
        ),
    ];

    let candidates = vec![
        Candidate::new(
            G01,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 1.0, None)],
        ),
        Candidate::new(
            G02,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 2.0, None)],
        ),
        Candidate::new(
            G03,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 3.0, None)],
        ),
        Candidate::new(
            G04,
            t,
            vec![Observation::pseudo_range(Carrier::L1, 4.0, None)],
        ),
    ];
}

#[test]
#[cfg(feature = "embed_ephem")]
fn cpp_matrix() {
    let cfg = Config::default()
        .with_modeling(Modeling::no_modeling())
        .with_navigation_method(Method::CPP);

    let almanac = crate::tests::almanac();
    let frame = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let coords_ecef_m = Vector3::new(
        REFERENCE_COORDS_ECEF_M.0,
        REFERENCE_COORDS_ECEF_M.1,
        REFERENCE_COORDS_ECEF_M.2,
    );

    let r_0 = (REFERENCE_COORDS_ECEF_M.0.powi(2)
        + REFERENCE_COORDS_ECEF_M.1.powi(2)
        + REFERENCE_COORDS_ECEF_M.2.powi(2))
    .sqrt();

    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let r0_orbit = Orbit::from_position(
        coords_ecef_m[0],
        coords_ecef_m[1],
        coords_ecef_m[2],
        t0_gpst,
        frame,
    );

    let initial_state = State::from_orbit(&r0_orbit)
        .unwrap();

    let mut candidates = vec![
        Candidate::new(
            G01,
            t0_gpst,
            vec![
                Observation::pseudo_range(Carrier::L1, 21401234.5, None),
                Observation::pseudo_range(Carrier::L2, 21401244.5, None),
            ],
        ),
        Candidate::new(
            G02,
            t0_gpst,
            vec![
                Observation::pseudo_range(Carrier::L1, 21401234.5, None),
                Observation::pseudo_range(Carrier::L2, 21401244.5, None),
            ],
        ),
        Candidate::new(
            G03,
            t0_gpst,
            vec![
                Observation::pseudo_range(Carrier::L1, 21401234.5, None),
                Observation::pseudo_range(Carrier::L2, 21401244.5, None),
            ],
        ),
        Candidate::new(
            G04,
            t0_gpst,
            vec![
                Observation::pseudo_range(Carrier::L1, 21401234.5, None),
                Observation::pseudo_range(Carrier::L2, 21401244.5, None),
            ],
        ),
    ];

    let null_bias = NullBias {};

    let null_time = NullTime {};
    let absolute_time = AbsoluteTime::new(null_time);

    let apriori = Apriori::from_ecef_m(coords_ecef_m, t0_gpst, frame);

    let sv_coords_m = vec![
        (15600.0, 7540.0, 20140.0),
        (18760.0, 2750.0, 18610.0),
        (17610.0, 14630.0, 13480.0),
        (19170.0, 610.0, 18390.0),
    ];

    for (nth, coords) in sv_coords_m.iter().enumerate() {
        let pos_vel_m = Vector6::new(coords.0, coords.1, coords.2, 0.0, 0.0, 0.0);
        let orbit = Orbit::from_cartesian_pos_vel(pos_vel_m / 1.0E3, t0_gpst, frame);
        candidates[nth].set_orbit(orbit);
        assert!(candidates[nth].orbit.is_some());
    }

    for cd in candidates.iter_mut() {
        cd.orbital_attitude_fixup(&almanac, r0_orbit).unwrap();
    }

    let mut nav = Navigation::new(
        &cfg,
        frame,
    );

    let r_i = vec![
        candidates[0].code_if_combination().unwrap().value,
        candidates[1].code_if_combination().unwrap().value,
        candidates[2].code_if_combination().unwrap().value,
        candidates[3].code_if_combination().unwrap().value,
    ];

    assert_eq!(
        r_i,
        vec![
            21401219.0427222,
            21401219.0427222,
            21401219.0427222,
            21401219.0427222
        ],
        "incorrect test values"
    );

    let rho = vec![
        ((REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[0].0).powi(2)
            + (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[0].1).powi(2)
            + (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[0].2).powi(2))
        .sqrt(),
        ((REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[1].0).powi(2)
            + (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[1].1).powi(2)
            + (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[1].2).powi(2))
        .sqrt(),
        ((REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[2].0).powi(2)
            + (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[2].1).powi(2)
            + (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[2].2).powi(2))
        .sqrt(),
        ((REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[3].0).powi(2)
            + (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[3].1).powi(2)
            + (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[3].2).powi(2))
        .sqrt(),
    ];

    // let (dx_m, dy_m, dz_m) = (
    //     (REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[i].0) / rho[i],
    //     (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[i].1) / rho[i],
    //     (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[i].2) / rho[i],
    // );

    // assert_eq!(nav.b[i], r_i[i] - rho[i], "b (noclock) test failed [{}]", i);

    nav.solving(t0_gpst, &initial_state, &candidates, 4, &null_bias, &absolute_time)
        .unwrap();

    // Clock definitions
    let mut cfg = Config::default().with_modeling(Modeling::no_modeling());

    cfg.modeling.sv_clock_bias = true;

    for (nth, cd) in candidates.iter_mut().enumerate() {
        let nanos = (nth + 100) as i128;
        let dt = Duration::from_total_nanoseconds(nanos);

        cd.set_clock_correction(ClockCorrection::without_relativistic_correction(dt));

        assert!(
            cd.clock_corr.is_some(),
            "failed to define a clock correction for {}@{}",
            cd.sv,
            cd.t
        );

        assert_eq!(
            cd.clock_corr.unwrap().duration.total_nanoseconds(),
            nanos,
            "invalid clock duration"
        );
    }

    let mut nav = Navigation::new(
        &cfg,
        frame,
    );

        // let (dx_m, dy_m, dz_m) = (
        //     (REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[i].0) / rho[i],
        //     (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[i].1) / rho[i],
        //     (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[i].2) / rho[i],
        // );

        // let dt_s = ((i + 100) as f64) * 1E-9;
        // assert_eq!(nav.b[i], r_i[i] - rho[i] + SPEED_OF_LIGHT_M_S * dt_s);

        nav.solving(t0_gpst, &initial_state, &candidates, 4, &null_bias, &absolute_time)
            .unwrap();

    cfg.modeling.relativistic_path_range = true;

    let mut nav = Navigation::new(
        &cfg,
        frame,
    );

    for i in 0..4 {
        let r_sat =
            (sv_coords_m[i].0.powi(2) + sv_coords_m[i].1.powi(2) + sv_coords_m[i].2.powi(2)).sqrt();
        let r_sat_0 = r_0 - r_sat;

        let dr = 2.0 * EARTH_GRAVITATION / SPEED_OF_LIGHT_M_S / SPEED_OF_LIGHT_M_S
            * ((r_sat + r_0 + r_sat_0) / (r_sat + r_0 - r_sat_0)).ln();

        // let (dx_m, dy_m, dz_m) = (
        //     (REFERENCE_COORDS_ECEF_M.0 - sv_coords_m[i].0) / (rho[i] + dr),
        //     (REFERENCE_COORDS_ECEF_M.1 - sv_coords_m[i].1) / (rho[i] + dr),
        //     (REFERENCE_COORDS_ECEF_M.2 - sv_coords_m[i].2) / (rho[i] + dr),
        // );

        // let dt_s = ((i + 100) as f64) * 1E-9;
        // let b_model = r_i[i] - (rho[i] + dr) + SPEED_OF_LIGHT_M_S * dt_s;

        // let err = (nav.b[i] - b_model).abs();
        // assert!(err < 1E-6);

        nav.solving(t0_gpst, &initial_state, &candidates, 4, &null_bias, &absolute_time)
            .unwrap();
    }
}
