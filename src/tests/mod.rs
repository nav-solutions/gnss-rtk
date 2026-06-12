use anise::{
    constants::frames::EARTH_J2000,
    prelude::{Almanac, Epoch, Frame, Orbit},
};

use std::str::FromStr;

use crate::navigation::apriori::Apriori;

pub mod ephemeris;
// pub mod fuzz;

mod bancroft;
mod candidate;
mod cpp;
mod data;
mod number;
mod phase_range;
mod pool;
mod ppp;
mod pseudo_range;
mod rtk;
mod spp;
mod time;

pub use data::*;

use log::LevelFilter;
use std::sync::Once;

pub const MAX_SPP_X_ERROR_M: f64 = 65.5;
pub const MAX_SPP_Y_ERROR_M: f64 = 39.3;
pub const MAX_SPP_Z_ERROR_M: f64 = 64.1;
pub const MAX_SPP_GDOP: f64 = 3.90;

pub const MAX_RTK_SPP_X_ERROR_M: f64 = 24.5;
pub const MAX_RTK_SPP_Y_ERROR_M: f64 = 11.5;
pub const MAX_RTK_SPP_Z_ERROR_M: f64 = 97.0;
pub const MAX_RTK_SPP_GDOP: f64 = 1.99;

pub const MAX_CPP_X_ERROR_M: f64 = 51.0;
pub const MAX_CPP_Y_ERROR_M: f64 = 38.0;
pub const MAX_CPP_Z_ERROR_M: f64 = 45.0;
pub const MAX_CPP_GDOP: f64 = 3.90;

pub const MAX_RTK_CPP_X_ERROR_M: f64 = 24.59;
pub const MAX_RTK_CPP_Y_ERROR_M: f64 = 11.5;
pub const MAX_RTK_CPP_Z_ERROR_M: f64 = 97.35;
pub const MAX_RTK_CPP_GDOP: f64 = 1.99;

pub const MAX_RTK_PPP_X_ERROR_M: f64 = 24.59;
pub const MAX_RTK_PPP_Y_ERROR_M: f64 = 11.5;
pub const MAX_RTK_PPP_Z_ERROR_M: f64 = 97.35;
pub const MAX_RTK_PPP_GDOP: f64 = 1.99;

pub const MAX_SURVEY_BANCROFT_X_ERROR_M: f64 = 35.0;
pub const MAX_SURVEY_BANCROFT_Y_ERROR_M: f64 = 110.0;
pub const MAX_SURVEY_BANCROFT_Z_ERROR_M: f64 = 30.0;

static INIT: Once = Once::new();

pub fn init_logger() {
    INIT.call_once(|| {
        env_logger::builder()
            .is_test(true)
            .filter_level(LevelFilter::Debug)
            .init();
    });
}

pub fn almanac() -> Almanac {
    Almanac::until_2035().unwrap_or_else(|e| panic!("Failed to build test Almanac: {e}"))
}

pub fn earth_frame() -> Frame {
    Almanac::until_2035()
        .unwrap_or_else(|e| panic!("Failed to build test Almanac: {e}"))
        .frame_from_uid(EARTH_J2000)
        .unwrap_or_else(|e| panic!("Failed to build test EARTH-J2000 frame: {e}"))
}

pub fn test_orbits() -> OrbitsData {
    let earth_frame = earth_frame();
    OrbitsData::new(earth_frame)
}

pub const ROVER_REFERENCE_COORDS_ECEF_M: (f64, f64, f64) =
    (3582105.2910, 532589.7313, 5232754.8054);

pub const BASE_REFERENCE_COORDS_ECEF_M: (f64, f64, f64) = (3628427.9118, 562059.0936, 5197872.2150);

pub fn reference_epoch() -> Epoch {
    Epoch::from_str("2020-06-25T00:00:00 GPST")
        .unwrap_or_else(|e| panic!("internal error: failed to build reference epoch: {e}"))
}

/// Builds reference [Apriori] position
pub fn rover_reference_apriori(t: Epoch) -> Apriori {
    let earth_frame = earth_frame();
    let ref_orbit = rover_reference_orbit(t, earth_frame);

    Apriori::from_orbit(&ref_orbit, earth_frame)
}

/// Builds reference [Apriori] position
pub fn base_reference_apriori(t: Epoch) -> Apriori {
    let earth_frame = earth_frame();
    let ref_orbit = base_reference_orbit(t, earth_frame);

    Apriori::from_orbit(&ref_orbit, earth_frame)
}

/// Builds reference [Apriori] position at reference [Epoch]
pub fn rover_reference_apriori_at_ref_epoch() -> Apriori {
    rover_reference_apriori(reference_epoch())
}

/// Builds reference [Apriori] position at reference [Epoch]
pub fn base_reference_apriori_at_ref_epoch() -> Apriori {
    base_reference_apriori(reference_epoch())
}

/// Expresses [ROVER_REFERENCE_COORDS_ECEF_M] as ANISE [Orbit] at any [Epoch].
pub fn rover_reference_orbit(t: Epoch, frame: Frame) -> Orbit {
    Orbit::from_position(
        ROVER_REFERENCE_COORDS_ECEF_M.0 / 1.0E3,
        ROVER_REFERENCE_COORDS_ECEF_M.1 / 1.0E3,
        ROVER_REFERENCE_COORDS_ECEF_M.2 / 1.0E3,
        t,
        frame,
    )
}

/// Expresses [BASE_REFERENCE_COORDS_ECEF_M] as ANISE [Orbit] at any [Epoch].
pub fn base_reference_orbit(t: Epoch, frame: Frame) -> Orbit {
    Orbit::from_position(
        BASE_REFERENCE_COORDS_ECEF_M.0 / 1.0E3,
        BASE_REFERENCE_COORDS_ECEF_M.1 / 1.0E3,
        BASE_REFERENCE_COORDS_ECEF_M.2 / 1.0E3,
        t,
        frame,
    )
}

/// Expresses [ROVER_REFERENCE_COORDS_ECEF_M] as ANISE [Orbit] at initial reference [Epoch].
pub fn rover_reference_orbit_at_ref_epoch(frame: Frame) -> Orbit {
    rover_reference_orbit(reference_epoch(), frame)
}

/// Expresses [BASE_REFERENCE_COORDS_ECEF_M] as ANISE [Orbit] at initial reference [Epoch].
pub fn base_reference_orbit_at_ref_epoch(frame: Frame) -> Orbit {
    base_reference_orbit(reference_epoch(), frame)
}

#[test]
fn verify_rover_reference_orbit() {
    let earth_frame = earth_frame();
    let orbit = rover_reference_orbit_at_ref_epoch(earth_frame);
    let apriori = rover_reference_apriori_at_ref_epoch();
    assert_eq!(orbit, apriori.to_orbit(), "invalid rover reference setup!");
}

#[test]
fn verify_rtkbase_reference_orbit() {
    let earth_frame = earth_frame();
    let orbit = base_reference_orbit_at_ref_epoch(earth_frame);
    let apriori = base_reference_apriori_at_ref_epoch();
    assert_eq!(
        orbit,
        apriori.to_orbit(),
        "invalid base station reference setup!"
    );
}
