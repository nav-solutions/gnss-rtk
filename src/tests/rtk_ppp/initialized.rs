use log::info;
use rstest::*;
use std::str::FromStr;

use crate::{
    navigation::apriori::Apriori,
    prelude::{Almanac, Config, Epoch, Frame, Method, Solver, UserParameters},
    tests::{
        CandidatesBuilder, OrbitsData, ROVER_REFERENCE_COORDS_ECEF_M, TestEnvironment,
        TestSpacebornBiases, ephemeris::NullEph, init_logger, time::NullTime,
    },
};

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
fn build_initial_apriori() -> Apriori {
    use crate::tests::rover_reference_apriori_at_ref_epoch;
    rover_reference_apriori_at_ref_epoch()
}

#[test]
fn static_rtk_ppp() {
    init_logger();

    let cfg = Config::default().with_navigation_method(Method::PPP);

    let default_params = UserParameters::default();

    let almanac = build_almanac();
    let earth_frame = build_earth_frame();

    let null_time = NullTime {};
    let null_eph = NullEph {};
    let environment = TestEnvironment::new();
    let space_biases = TestSpacebornBiases::build();

    let orbits_data = OrbitsData::new(earth_frame);

    let rtk_base = CandidatesBuilder::build_rtk_base();

    let mut solver = Solver::new(
        almanac,
        earth_frame,
        cfg,
        null_eph.into(),
        orbits_data.into(),
        space_biases.into(),
        environment.into(),
        null_time,
        Some(ROVER_REFERENCE_COORDS_ECEF_M),
    );

    for (nth, epoch_str) in [
        "2020-06-25T00:00:00 GPST",
        "2020-06-25T00:15:00 GPST",
        // "2020-06-25T00:30:00 GPST", //TODO RTK test data
        // "2020-06-25T00:45:00 GPST", // TODO RTK test data
        // "2020-06-25T01:00:00 GPST", // TODO RTK test data
    ]
    .iter()
    .enumerate()
    {
        let t_gpst = Epoch::from_str(epoch_str).unwrap();

        let candidates = CandidatesBuilder::build_rover_at(t_gpst);

        assert!(
            !candidates.is_empty(),
            "{epoch_str} - no measurements to propose"
        );

        let status = solver.rtk(t_gpst, default_params, &candidates, &rtk_base);

        match status {
            Err(e) => panic!("Static RTK-PPP process failed with invalid error: {e}"),
            Ok(pvt) => {
                info!("Solution #{} {:#?}", nth + 1, pvt);

                let (pos_x_m, pos_y_m, pos_z_m) = pvt.pos_m;
                let (expected_x_m, expected_y_m, expected_z_m) = ROVER_REFERENCE_COORDS_ECEF_M;

                let (err_x_m, err_y_m, err_z_m) = (
                    (pos_x_m - expected_x_m).abs(),
                    (pos_y_m - expected_y_m).abs(),
                    (pos_z_m - expected_z_m).abs(),
                );

                assert!(
                    err_x_m < 100.0,
                    "epoch={epoch_str} - x error={err_x_m}m too large"
                );

                assert!(
                    err_y_m < 100.0,
                    "epoch={epoch_str} - y error={err_y_m}m too large"
                );

                assert!(
                    err_z_m < 100.0,
                    "epoch={epoch_str} - z error={err_z_m}m too large"
                );

                info!(
                    "{epoch_str} (static) rtk-ppp (with preset) error: x={err_x_m}m y={err_y_m}m z={err_z_m}"
                );
            },
        }
    }
}
