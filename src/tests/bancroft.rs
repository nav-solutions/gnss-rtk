use crate::{
    bancroft::Bancroft,
    prelude::{Epoch, Frame, OrbitSource, SpacebornBias},
    tests::{
        CandidatesBuilder, MAX_SURVEY_BANCROFT_X_ERROR_M, MAX_SURVEY_BANCROFT_Y_ERROR_M,
        MAX_SURVEY_BANCROFT_Z_ERROR_M, OrbitsData, ROVER_REFERENCE_COORDS_ECEF_M,
        TestSpacebornBiases, init_logger,
    },
};

use log::info;
use rstest::*;

use std::str::FromStr;

#[fixture]
fn build_earth_frame() -> Frame {
    use crate::tests::earth_frame;
    earth_frame()
}

#[test]
fn bancroft_gpst() {
    init_logger();

    let earth_frame = build_earth_frame();
    let orbits_data = OrbitsData::new(earth_frame);
    let space_biases = TestSpacebornBiases::build();

    for t_str in [
        "2020-06-25T00:00:00 GPST",
        "2020-06-25T00:15:00 GPST",
        "2020-06-25T00:30:00 GPST",
        "2020-06-25T00:45:00 GPST",
        "2020-06-25T01:00:00 GPST",
    ] {
        let epoch = Epoch::from_str(t_str).unwrap_or_else(|e| {
            panic!("invalid test epoch {t_str}: {e}");
        });

        let mut pool = CandidatesBuilder::build_rover_at(epoch);

        for cd in pool.iter_mut() {
            let state = orbits_data
                .state_at(cd.epoch, cd.sv, earth_frame)
                .unwrap_or_else(|| {
                    panic!("{}({}) - orbital state not defined", cd.epoch, cd.sv);
                });

            cd.orbit = Some(state);

            let rtm = cd.to_partial_bias_runtime();

            cd.clock_corr = space_biases.clock_bias(&rtm);
            cd.tgd = space_biases.group_delay(&rtm);
        }

        let solver = Bancroft::new(&pool)
            .unwrap_or_else(|e| panic!("failed to create Bancroft solver: {e}"));

        let output = solver
            .resolve()
            .unwrap_or_else(|e| panic!("Bancroft solver failure: {e}"));

        let x_err = (output[0] - ROVER_REFERENCE_COORDS_ECEF_M.0).abs();
        let y_err = (output[1] - ROVER_REFERENCE_COORDS_ECEF_M.1).abs();
        let z_err = (output[2] - ROVER_REFERENCE_COORDS_ECEF_M.2).abs();

        assert!(
            x_err < MAX_SURVEY_BANCROFT_X_ERROR_M,
            "{t_str} bancroft x-error too large: {x_err}"
        );
        assert!(
            y_err < MAX_SURVEY_BANCROFT_Y_ERROR_M,
            "{t_str} bancroft y-error too large: {y_err}"
        );
        assert!(
            z_err < MAX_SURVEY_BANCROFT_Z_ERROR_M,
            "{t_str} bancroft z-error too large: {z_err}"
        );

        info!("{t_str} bancroft errors x={x_err} y={y_err} z={z_err}");
    }
}
