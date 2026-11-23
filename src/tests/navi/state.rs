use crate::{
    navigation::State,
    prelude::{Almanac, Epoch, Orbit, Vector3, EARTH_J2000, SPEED_OF_LIGHT_M_S},
    tests::data::gps::J2020_06_25_GPS_EPOCHS,
};

use anise::math::{Vector4, Vector6};
use std::str::FromStr;

#[test]
fn state_initialization() {
    let almanac = crate::tests::almanac();
    let frame = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let t0 = Epoch::from_str(J2020_06_25_GPS_EPOCHS[0]).unwrap();

    let coords = Vector3::new(1.0, 2.0, 3.0);

    let mut state = State::from_ecef_m(coords, t0, frame).unwrap();

    assert_eq!(state.pos_m, (1.0, 2.0, 3.0));
    assert_eq!(state.dt.total_nanoseconds(), 0);

    let dx = Vector4::new(1.0, 2.0, 3.0, 4.0);
    state.update(dx);

    assert_eq!(state.pos_m, (2.0, 4.0, 6.0));

    let dt_nanos = (4.0 / SPEED_OF_LIGHT_M_S * 1E9).round() as i128;
    assert_eq!(state.dt.total_nanoseconds(), dt_nanos);

    let pos_vel_km = Vector6::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
    let orbit = Orbit::from_cartesian_pos_vel(pos_vel_km, t0, frame);

    let state = State::from_orbit(&orbit).unwrap();

    assert_eq!(state.pos_m, (1000.0, 2000.0, 3000.0));
    assert_eq!(state.dt.total_nanoseconds(), 0);
}
