use log::info;
use rstest::*;
use std::str::FromStr;

use crate::{
    navigation::apriori::Apriori,
    prelude::{
        Almanac, ClockProfile, Config, Epoch, Frame, Method, PVTSolutionType, Solver,
        UserParameters, UserProfile,
    },
    tests::{
        CandidatesBuilder, MAX_SPP_GDOP, MAX_SPP_X_ERROR_M, MAX_SPP_Y_ERROR_M, MAX_SPP_Z_ERROR_M,
        OrbitsData, ROVER_REFERENCE_COORDS_ECEF_M, TestEnvironment, TestSpacebornBiases,
        ephemeris::NullEph, init_logger, time::NullTime,
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
fn static_spp() {
    init_logger();

    let cfg = Config::default().with_navigation_method(Method::SPP);

    let default_params = UserParameters::new(UserProfile::Static, ClockProfile::Quartz);

    let almanac = build_almanac();
    let earth_frame = build_earth_frame();

    let null_time = NullTime {};
    let null_eph = NullEph {};
    let environment = TestEnvironment::new();
    let space_biases = TestSpacebornBiases::build();

    let orbits_data = OrbitsData::new(earth_frame);

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
        "2020-06-25T00:30:00 GPST",
        "2020-06-25T00:45:00 GPST",
        //"2020-06-25T01:00:00 GPST",
    ]
    .iter()
    .enumerate()
    {
        let t_gpst = Epoch::from_str(epoch_str).unwrap();

        let candidates = CandidatesBuilder::build_rover_at(t_gpst);

        assert!(
            !candidates.is_empty(),
            "no measurements to propose at \"{epoch_str}\""
        );

        let status = solver.ppp(t_gpst, default_params, &candidates);

        match status {
            Err(e) => panic!("Static SPP process failed with invalid error: {e}"),
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
                    err_x_m < MAX_SPP_X_ERROR_M,
                    "epoch={epoch_str} - x error={err_x_m:.3}m too large"
                );

                assert!(
                    err_y_m < MAX_SPP_Y_ERROR_M,
                    "epoch={epoch_str} - y error={err_y_m:.3}m too large"
                );

                assert!(
                    err_z_m < MAX_SPP_Z_ERROR_M,
                    "epoch={epoch_str} - z error={err_z_m:.3}m too large"
                );

                assert_eq!(pvt.solution_type, PVTSolutionType::PPP);

                assert!(
                    pvt.gdop < MAX_SPP_GDOP,
                    "{epoch_str} (static) spp survey GDOP too large!"
                );

                info!(
                    "{} (static) spp (with preset) error: x={:.3}m y={:.3}m z={:.3}m, GDOP={} TDOP={}",
                    epoch_str, err_x_m, err_y_m, err_z_m, pvt.gdop, pvt.tdop,
                );
            },
        }
    }
}
