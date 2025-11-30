use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    navigation::{apriori::Apriori, state::State},
    pool::Pool,
    prelude::{Almanac, Carrier, Config, Epoch, Frame, Method, TimeScale},
    tests::{
        E01, E03, E05, OrbitsData, TestEnvironment, TestSpacebornBiases, data::CandidatesBuilder,
        ephemeris::NullEph, init_logger, time::NullTime,
    },
};

use rstest::*;

use std::rc::Rc;

use hifitime::Duration;

use std::str::FromStr;

#[fixture]
fn build_almanac() -> Almanac {
    use crate::tests::almanac;
    almanac()
}

#[fixture]
fn build_earth_frame() -> Frame {
    use crate::tests::earth_frame;
    earth_frame()
}

#[fixture]
fn build_orbit_source() -> OrbitsData {
    use crate::tests::test_orbits;
    test_orbits()
}

#[fixture]
fn build_rover_apriori() -> Apriori {
    use crate::tests::rover_reference_apriori_at_ref_epoch;
    rover_reference_apriori_at_ref_epoch()
}

#[fixture]
fn build_base_apriori() -> Apriori {
    use crate::tests::base_reference_apriori_at_ref_epoch;
    base_reference_apriori_at_ref_epoch()
}

#[test]
fn ppp_pool_fit() {
    init_logger();

    let almanac = build_almanac();
    let default_cfg = Config::default();

    let null_eph = NullEph {};
    let null_time = NullTime {};

    let earth_frame = build_earth_frame();
    let orbits_data = Rc::new(build_orbit_source());
    let environment = TestEnvironment::new();
    let space_biases = TestSpacebornBiases::build();

    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let candidates = CandidatesBuilder::build_rover_at(t0_gpst);

    let mut pool = Pool::allocate(
        almanac,
        default_cfg,
        earth_frame,
        null_eph.into(),
        orbits_data,
        environment.into(),
        space_biases.into(),
    );

    pool.new_epoch(&candidates);

    assert_eq!(
        pool.candidates().len(),
        candidates.len(),
        "builder dropped some data!"
    );

    pool.pre_fit("rover", &null_time);
    pool.orbital_states_fit("rover");

    assert_eq!(
        pool.candidates().len(),
        candidates.len(),
        "pool fit dropped some data!"
    );
}

#[test]
fn ppp_gst_pool_fit() {
    init_logger();

    let almanac = build_almanac();
    let earth_frame = build_earth_frame();

    let mut default_cfg = Config::default();
    default_cfg.timescale = TimeScale::GST;

    let null_eph = NullEph {};
    let null_time = NullTime {};

    let orbits_data = Rc::new(build_orbit_source());
    let environment = TestEnvironment::new();
    let space_biases = TestSpacebornBiases::build();

    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let candidates = CandidatesBuilder::build_rover_at(t0_gpst);

    let mut pool = Pool::allocate(
        almanac,
        default_cfg,
        earth_frame,
        null_eph.into(),
        orbits_data,
        environment.into(),
        space_biases.into(),
    );

    pool.new_epoch(&candidates);
    pool.pre_fit("rover", &null_time);

    assert_eq!(
        pool.candidates().len(),
        candidates.len(),
        "pool prefit dropped some data!"
    );

    for cd in pool.candidates().iter() {
        assert_eq!(
            cd.epoch.time_scale,
            TimeScale::GST,
            "{}({}) - not correctly transposed to GST",
            cd.epoch,
            cd.sv
        );
    }

    pool.orbital_states_fit("rover");

    for cd in pool.candidates().iter() {
        assert_eq!(
            cd.epoch.time_scale,
            TimeScale::GST,
            "{}({}) - not correctly transposed to GST",
            cd.epoch,
            cd.sv
        );

        assert_eq!(
            cd.tx_epoch.time_scale,
            TimeScale::GST,
            "{}({}) - not correctly transposed to GST",
            cd.tx_epoch,
            cd.sv
        );

        assert!(
            cd.epoch > cd.tx_epoch,
            "{}({}) - physical non sense (after transposition to GST)",
            cd.epoch,
            cd.sv
        );

        let dt = cd.epoch - cd.tx_epoch;

        assert!(
            dt < Duration::from_milliseconds(100.0),
            "{}({}) - not compatible with MEO/LEO Earth nav (after transposition to GST)",
            cd.epoch,
            cd.sv
        );
    }
}

#[test]
fn rtk_spp_pool_fit() {
    init_logger();

    let almanac = build_almanac();

    let default_cfg = Config::default().with_navigation_method(Method::SPP);

    let null_time = NullTime {};
    let null_eph = Rc::new(NullEph {});

    let earth_frame = build_earth_frame();
    let orbits_data = Rc::new(build_orbit_source());
    let environment = Rc::new(TestEnvironment::new());
    let space_biases = Rc::new(TestSpacebornBiases::build());

    const NUM_SV: usize = 8;
    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let mut rover = Pool::allocate(
        almanac.clone(),
        default_cfg.clone(),
        earth_frame,
        null_eph.clone(),
        orbits_data.clone(),
        environment.clone(),
        space_biases.clone(),
    );

    let candidates = CandidatesBuilder::build_rover_at(t0_gpst);

    rover.new_epoch(&candidates);

    let mut base = Pool::allocate(
        almanac,
        default_cfg,
        earth_frame,
        null_eph,
        orbits_data,
        environment,
        space_biases,
    );

    let candidates = CandidatesBuilder::build_base_at(t0_gpst);

    base.new_epoch(&candidates);

    let apriori = build_rover_apriori();

    let rover_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build rover initial state: {e}");
    });

    let apriori = build_base_apriori();

    let base_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build base initial state: {e}");
    });

    base.pre_fit("base", &null_time);
    base.orbital_states_fit("base");

    base.post_fit("base", &base_state).unwrap_or_else(|e| {
        panic!("base station post-fit failed with {e}");
    });

    rover.pre_fit("rover", &null_time);
    rover.orbital_states_fit("rover");

    rover.post_fit("rover", &rover_state).unwrap_or_else(|e| {
        panic!("rover post-fit failed with {e}");
    });

    let dbl_diff = rover.rtk_post_fit(&mut base).unwrap_or_else(|e| {
        panic!("rtk post-fit failed with {e}");
    });

    let expected_sat = NUM_SV - 2; // pivot+E13 below elevation mask

    assert!(
        dbl_diff.inner.len() == expected_sat,
        "Did not form correct DDs",
    );

    let mut e01_passed = false;
    let mut e03_passed = false;

    for (sat, dd) in dbl_diff.inner.iter() {
        assert!(*sat != E05, "did not remove pivot sat!");

        if *sat == E01 {
            let (carrier, code) = dd.code.expect("E01 missing DD(code)");

            assert_eq!(carrier, Carrier::L1);

            assert_eq!(
                code,
                27616185.992 - 23730317.923 - (27506424.743 - 23595077.027)
            );

            e01_passed = true;
        } else if *sat == E03 {
            let (carrier, code) = dd.code.expect("E03 missing DD(code)");

            assert_eq!(carrier, Carrier::L1);

            assert_eq!(
                code,
                27055946.391 - 23730317.923 - (26952639.751 - 23595077.027)
            );

            e03_passed = true;
        }
    }

    assert!(e01_passed, "E01 test failed");
    assert!(e03_passed, "E03 test failed");
}

#[test]
fn rtk_cpp_pool_fit() {
    init_logger();

    let almanac = build_almanac();

    let default_cfg = Config::default().with_navigation_method(Method::CPP);

    let null_time = NullTime {};
    let null_eph = Rc::new(NullEph {});

    let earth_frame = build_earth_frame();
    let orbits_data = Rc::new(build_orbit_source());
    let environment = Rc::new(TestEnvironment::new());
    let space_biases = Rc::new(TestSpacebornBiases::build());

    const NUM_SV: usize = 8;
    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let mut rover = Pool::allocate(
        almanac.clone(),
        default_cfg.clone(),
        earth_frame,
        null_eph.clone(),
        orbits_data.clone(),
        environment.clone(),
        space_biases.clone(),
    );

    let candidates = CandidatesBuilder::build_rover_at(t0_gpst);

    rover.new_epoch(&candidates);

    let mut base = Pool::allocate(
        almanac,
        default_cfg,
        earth_frame,
        null_eph,
        orbits_data,
        environment,
        space_biases,
    );

    let candidates = CandidatesBuilder::build_base_at(t0_gpst);

    base.new_epoch(&candidates);

    let apriori = build_rover_apriori();

    let rover_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build rover initial state: {e}");
    });

    let apriori = build_base_apriori();

    let base_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build base initial state: {e}");
    });

    base.pre_fit("base", &null_time);
    base.orbital_states_fit("base");

    base.post_fit("base", &base_state).unwrap_or_else(|e| {
        panic!("base station post-fit failed with {e}");
    });

    rover.pre_fit("rover", &null_time);
    rover.orbital_states_fit("rover");

    rover.post_fit("rover", &rover_state).unwrap_or_else(|e| {
        panic!("rover post-fit failed with {e}");
    });

    let dbl_diff = rover.rtk_post_fit(&mut base).unwrap_or_else(|e| {
        panic!("rtk post-fit failed with {e}");
    });

    let expected_sat = NUM_SV - 2; // pivot+E13 below elevation mask

    assert!(
        dbl_diff.inner.len() == expected_sat,
        "did not form correct DD observations"
    );

    let mut e01_passed = false;
    let mut e03_passed = false;

    let (f1, f2) = (Carrier::L1.frequency_hz(), Carrier::E5b.frequency_hz());
    let (f1pow, f2pow) = (f1.powi(2), f2.powi(2));

    let pc_e01_rover = (f1pow * 27616185.992 - f2pow * 27616184.819) / (f1pow - f2pow);
    let pc_e03_rover = (f1pow * 27055946.391 - f2pow * 27055945.532) / (f1pow - f2pow);
    let pc_e05_rover = (f1pow * 23730317.923 - f2pow * 23730316.788) / (f1pow - f2pow);

    let pc_e01_base = (f1pow * 27506424.743 - f2pow * 27506425.902) / (f1pow - f2pow);
    let pc_e03_base = (f1pow * 26952639.751 - f2pow * 26952641.150) / (f1pow - f2pow);
    let pc_e05_base = (f1pow * 23595077.027 - f2pow * 23595078.180) / (f1pow - f2pow);

    for (sat, dd) in dbl_diff.inner.iter() {
        assert!(*sat != E05, "did not remove pivot sat!");

        if *sat == E01 {
            let (carrier, code) = dd.code_if.expect("E01 missing DD(code_if)");

            assert_eq!(carrier, Carrier::L1);

            assert_eq!(
                code,
                pc_e01_rover - pc_e05_rover - (pc_e01_base - pc_e05_base)
            );

            e01_passed = true;
        } else if *sat == E03 {
            let (carrier, code) = dd.code_if.expect("E03 missing DD(code_if)");

            assert_eq!(carrier, Carrier::L1);
            assert_eq!(
                code,
                pc_e03_rover - pc_e05_rover - (pc_e03_base - pc_e05_base)
            );

            e03_passed = true;
        }
    }

    assert!(e01_passed, "E01 test failed");
    assert!(e03_passed, "E03 test failed");
}

#[test]
fn rtk_ppp_pool_fit() {
    init_logger();

    let almanac = build_almanac();

    let default_cfg = Config::default().with_navigation_method(Method::PPP);

    let null_time = NullTime {};
    let null_eph = Rc::new(NullEph {});

    let earth_frame = build_earth_frame();
    let orbits_data = Rc::new(build_orbit_source());
    let environment = Rc::new(TestEnvironment::new());
    let space_biases = Rc::new(TestSpacebornBiases::build());

    const NUM_SV: usize = 8;
    let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

    let mut rover = Pool::allocate(
        almanac.clone(),
        default_cfg.clone(),
        earth_frame,
        null_eph.clone(),
        orbits_data.clone(),
        environment.clone(),
        space_biases.clone(),
    );

    let candidates = CandidatesBuilder::build_rover_at(t0_gpst);

    rover.new_epoch(&candidates);

    let mut base = Pool::allocate(
        almanac,
        default_cfg,
        earth_frame,
        null_eph,
        orbits_data,
        environment,
        space_biases,
    );

    let candidates = CandidatesBuilder::build_base_at(t0_gpst);

    base.new_epoch(&candidates);

    let apriori = build_rover_apriori();

    let rover_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build rover initial state: {e}");
    });

    let apriori = build_base_apriori();

    let base_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
        panic!("Failed to build base initial state: {e}");
    });

    base.pre_fit("base", &null_time);
    base.orbital_states_fit("base");

    base.post_fit("base", &base_state).unwrap_or_else(|e| {
        panic!("base station post-fit failed with {e}");
    });

    rover.pre_fit("rover", &null_time);
    rover.orbital_states_fit("rover");

    rover.post_fit("rover", &rover_state).unwrap_or_else(|e| {
        panic!("rover post-fit failed with {e}");
    });

    let dbl_diff = rover.rtk_post_fit(&mut base).unwrap_or_else(|e| {
        panic!("rtk post-fit failed with {e}");
    });

    let expected_sat = NUM_SV - 2; // pivot+E13 below elevation mask

    assert!(
        dbl_diff.inner.len() == expected_sat,
        "did not form correct DD observations"
    );

    let mut e01_passed = false;
    let mut e03_passed = false;

    let p11_rover = 27616185.992;
    let p12_rover = 27616184.819;
    let p31_rover = 27055946.391;
    let p32_rover = 27055945.532;
    let p51_rover = 23730317.923;
    let p52_rover = 23730316.788;

    let p11_base = 27506424.743;
    let p12_base = 27506425.902;
    let p31_base = 26952639.751;
    let p32_base = 26952641.150;
    let p51_base = 23595077.027;
    let p52_base = 23595078.180;

    let l11_rover = 145124050.106;
    let l12_rover = 108371872.760;
    let l31_rover = 142179967.778;
    let l32_rover = 106173364.686;
    let l51_rover = 124703702.220;
    let l52_rover = 93122915.921;

    let l11_base = 144547262.730;
    let l12_base = 107941157.967;
    let l31_base = 141637100.045;
    let l32_base = 105767978.933;
    let l51_base = 123992999.980;
    let l52_base = 92592203.174;

    let (f1, f2) = (Carrier::L1.frequency_hz(), Carrier::E5b.frequency_hz());
    let (f1pow, f2pow) = (f1.powi(2), f2.powi(2));

    let f_if = f1 * f2 / (f1pow + f2pow).sqrt();
    let lambda_if = SPEED_OF_LIGHT_M_S / f_if;

    let (lambda_n, lambda_w) = (
        SPEED_OF_LIGHT_M_S / (f1 + f2),
        SPEED_OF_LIGHT_M_S / (f1 - f2),
    );

    let pc_e01_rover = (f1pow * p11_rover - f2pow * p12_rover) / (f1pow - f2pow);
    let pc_e03_rover = (f1pow * p31_rover - f2pow * p32_rover) / (f1pow - f2pow);
    let pc_e05_rover = (f1pow * p51_rover - f2pow * p52_rover) / (f1pow - f2pow);
    let lc_e01_rover = (f1pow * l11_rover - f2pow * l12_rover) / (f1pow - f2pow);
    let lc_e03_rover = (f1pow * l31_rover - f2pow * l32_rover) / (f1pow - f2pow);
    let lc_e05_rover = (f1pow * l31_rover - f2pow * l52_rover) / (f1pow - f2pow);

    let pc_e01_base = (f1pow * p11_base - f2pow * p12_base) / (f1pow - f2pow);
    let pc_e03_base = (f1pow * p31_base - f2pow * p32_base) / (f1pow - f2pow);
    let pc_e05_base = (f1pow * p51_base - f2pow * p52_base) / (f1pow - f2pow);
    let lc_e01_base = (f1pow * l11_base - f2pow * l12_base) / (f1pow - f2pow);
    let lc_e03_base = (f1pow * l31_base - f2pow * l32_base) / (f1pow - f2pow);
    let lc_e05_base = (f1pow * l31_base - f2pow * l52_base) / (f1pow - f2pow);

    for (sat, dd) in dbl_diff.inner.iter() {
        assert!(*sat != E05, "pivot sat should have been removed.");

        if *sat == E01 {
            assert_eq!(
                dd.code,
                Some((Carrier::L1, p11_rover - p51_rover - p11_base + p51_base)),
                "{}(E01): invalid DD(C1)",
                t0_gpst,
            );

            assert_eq!(
                dd.code_j,
                Some((Carrier::E5b, p12_rover - p52_rover - p12_base + p52_base)),
                "{}(E01): invalid DD(C2)",
                t0_gpst,
            );

            let pc_e01_rover = (f1pow * p11_rover - f2pow * p12_rover) / (f1pow - f2pow);
            let pc_e05_rover = (f1pow * p51_rover - f2pow * p52_rover) / (f1pow - f2pow);
            let pc_e01_base = (f1pow * p11_base - f2pow * p12_base) / (f1pow - f2pow);
            let pc_e05_base = (f1pow * p51_base - f2pow * p52_base) / (f1pow - f2pow);

            assert_eq!(
                dd.code_if,
                Some((
                    Carrier::L1,
                    pc_e01_rover - pc_e05_rover - pc_e01_base + pc_e05_base
                )),
                "{}(E01): invalid DD(C2)",
                t0_gpst,
            );

            assert_eq!(
                dd.phase,
                Some((Carrier::L1, l11_rover - l51_rover - l11_base + l51_base)),
                "{}(E01): invalid DD(L1)",
                t0_gpst,
            );

            assert_eq!(
                dd.phase_j,
                Some((Carrier::E5b, l12_rover - l52_rover - l12_base + l52_base)),
                "{}(E01): invalid DD(L2)",
                t0_gpst,
            );

            let lc_e01_rover = (f1pow * l11_rover - f2pow * l12_rover) / (f1pow - f2pow);
            let lc_e05_rover = (f1pow * l51_rover - f2pow * l52_rover) / (f1pow - f2pow);
            let lc_e01_base = (f1pow * l11_base - f2pow * l12_base) / (f1pow - f2pow);
            let lc_e05_base = (f1pow * l51_base - f2pow * l52_base) / (f1pow - f2pow);

            assert_eq!(
                dd.phase_if,
                Some((
                    Carrier::L1,
                    lambda_if,
                    lc_e01_rover - lc_e05_rover - lc_e01_base + lc_e05_base
                )),
                "{}(E01): invalid DD(Lif)",
                t0_gpst,
            );

            let lw_e01_rover = (f1 * l11_rover - f2 * l12_rover) / (f1 - f2);
            let lw_e05_rover = (f1 * l51_rover - f2 * l52_rover) / (f1 - f2);
            let lw_e01_base = (f1 * l11_base - f2 * l12_base) / (f1 - f2);
            let lw_e05_base = (f1 * l51_base - f2 * l52_base) / (f1 - f2);

            assert_eq!(
                dd.lw,
                Some((
                    Carrier::L1,
                    lambda_w,
                    lw_e01_rover - lw_e05_rover - lw_e01_base + lw_e05_base
                )),
                "{}(E01): invalid DD(Lw)",
                t0_gpst,
            );

            let cn_e01_rover = (f1 * p11_rover + f2 * p12_rover) / (f1 + f2);
            let cn_e05_rover = (f1 * p51_rover + f2 * p52_rover) / (f1 + f2);
            let cn_e01_base = (f1 * p11_base + f2 * p12_base) / (f1 + f2);
            let cn_e05_base = (f1 * p51_base + f2 * p52_base) / (f1 + f2);

            assert_eq!(
                dd.cn,
                Some((
                    Carrier::L1,
                    lambda_n,
                    cn_e01_rover - cn_e05_rover - cn_e01_base + cn_e05_base
                )),
                "{}(E01): invalid DD(Cn)",
                t0_gpst,
            );

            e01_passed = true;
        } else if *sat == E03 {
            assert_eq!(
                dd.code,
                Some((Carrier::L1, p31_rover - p51_rover - p31_base + p51_base)),
                "{}(E03): invalid DD(C1)",
                t0_gpst,
            );

            assert_eq!(
                dd.code_j,
                Some((Carrier::E5b, p32_rover - p52_rover - p32_base + p52_base)),
                "{}(E03): invalid DD(C2)",
                t0_gpst,
            );

            let pc_e03_rover = (f1pow * p31_rover - f2pow * p32_rover) / (f1pow - f2pow);
            let pc_e05_rover = (f1pow * p51_rover - f2pow * p52_rover) / (f1pow - f2pow);
            let pc_e03_base = (f1pow * p31_base - f2pow * p32_base) / (f1pow - f2pow);
            let pc_e05_base = (f1pow * p51_base - f2pow * p52_base) / (f1pow - f2pow);

            assert_eq!(
                dd.code_if,
                Some((
                    Carrier::L1,
                    pc_e03_rover - pc_e05_rover - pc_e03_base + pc_e05_base
                )),
                "{}(E03): invalid DD(C2)",
                t0_gpst,
            );

            assert_eq!(
                dd.phase,
                Some((Carrier::L1, l31_rover - l51_rover - l31_base + l51_base)),
                "{}(E03): invalid DD(L1)",
                t0_gpst,
            );

            assert_eq!(
                dd.phase_j,
                Some((Carrier::E5b, l32_rover - l52_rover - l32_base + l52_base)),
                "{}(E03): invalid DD(L2)",
                t0_gpst,
            );

            let lc_e03_rover = (f1pow * l31_rover - f2pow * l32_rover) / (f1pow - f2pow);
            let lc_e05_rover = (f1pow * l51_rover - f2pow * l52_rover) / (f1pow - f2pow);
            let lc_e01_base = (f1pow * l31_base - f2pow * l32_base) / (f1pow - f2pow);
            let lc_e05_base = (f1pow * l51_base - f2pow * l52_base) / (f1pow - f2pow);

            assert_eq!(
                dd.phase_if,
                Some((
                    Carrier::L1,
                    lambda_if,
                    lc_e03_rover - lc_e05_rover - lc_e03_base + lc_e05_base
                )),
                "{}(E03): invalid DD(Lif)",
                t0_gpst,
            );

            let lw_e03_rover = (f1 * l31_rover - f2 * l32_rover) / (f1 - f2);
            let lw_e05_rover = (f1 * l51_rover - f2 * l52_rover) / (f1 - f2);
            let lw_e03_base = (f1 * l31_base - f2 * l32_base) / (f1 - f2);
            let lw_e05_base = (f1 * l51_base - f2 * l52_base) / (f1 - f2);

            assert_eq!(
                dd.lw,
                Some((
                    Carrier::L1,
                    lambda_w,
                    lw_e03_rover - lw_e05_rover - lw_e03_base + lw_e05_base
                )),
                "{}(E03): invalid DD(Lw)",
                t0_gpst,
            );

            let cn_e03_rover = (f1 * p31_rover + f2 * p32_rover) / (f1 + f2);
            let cn_e05_rover = (f1 * p51_rover + f2 * p52_rover) / (f1 + f2);
            let cn_e03_base = (f1 * p31_base + f2 * p32_base) / (f1 + f2);
            let cn_e05_base = (f1 * p51_base + f2 * p52_base) / (f1 + f2);

            assert_eq!(
                dd.cn,
                Some((
                    Carrier::L1,
                    lambda_n,
                    cn_e03_rover - cn_e05_rover - cn_e03_base + cn_e05_base
                )),
                "{}(E03): invalid DD(Cn)",
                t0_gpst,
            );

            e03_passed = true;
        }
    }

    assert!(e01_passed, "E01 test failed");
    assert!(e03_passed, "E03 test failed");
}
