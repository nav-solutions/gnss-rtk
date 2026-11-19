use anise::{
    almanac::metaload::{MetaAlmanac, MetaFile},
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

pub mod sim;

pub use data::*;
mod almanac;

pub use almanac::{almanac, earth_frame};

use log::LevelFilter;
use std::{
    path::{Path, PathBuf},
    sync::{Once, OnceLock},
};

// Performance budgets of the frozen data tests (ESBC00DNK rover, MOJN00DNK base,
// 2020-06-25, Galileo, GRG SP3 orbits + 30s clocks), absolute ECEF errors per axis.
//
// Budgets are about 1.3 to 1.5 x the worst error observed on the frozen epochs:
// the processing is deterministic, so any degradation is caught. Tighten them
// when the solver improves (each budget documents what limits it today).
//
// SPP: L1 only, the test environment does not model the ionosphere (~30 m vertical).
pub const MAX_SPP_X_ERROR_M: f64 = 15.0;
pub const MAX_SPP_Y_ERROR_M: f64 = 5.0;
pub const MAX_SPP_Z_ERROR_M: f64 = 35.0;
pub const MAX_SPP_GDOP: f64 = 3.90;

// RTK-SPP: 35 km baseline, single frequency double differences.
pub const MAX_RTK_SPP_X_ERROR_M: f64 = 0.6;
pub const MAX_RTK_SPP_Y_ERROR_M: f64 = 0.5;
pub const MAX_RTK_SPP_Z_ERROR_M: f64 = 0.45;
pub const MAX_RTK_SPP_GDOP: f64 = 1.99;

// CPP: iono free code. Limited by code noise and unmodeled antenna offsets
// (reference coordinates are the marker, the ARP is 21.6 cm above).
pub const MAX_CPP_X_ERROR_M: f64 = 0.8;
pub const MAX_CPP_Y_ERROR_M: f64 = 0.8;
pub const MAX_CPP_Z_ERROR_M: f64 = 2.2;
pub const MAX_CPP_GDOP: f64 = 3.90;

// RTK-CPP: also used by the RTK -> CPP context switch test (absolute epochs).
pub const MAX_RTK_CPP_X_ERROR_M: f64 = 0.8;
pub const MAX_RTK_CPP_Y_ERROR_M: f64 = 0.6;
pub const MAX_RTK_CPP_Z_ERROR_M: f64 = 1.2;
pub const MAX_RTK_CPP_GDOP: f64 = 1.99;

pub const MAX_RTK_PPP_X_ERROR_M: f64 = 0.25;
pub const MAX_RTK_PPP_Y_ERROR_M: f64 = 0.65;
pub const MAX_RTK_PPP_Z_ERROR_M: f64 = 0.2;
pub const MAX_RTK_PPP_GDOP: f64 = 1.99;

// Bancroft: raw L1 pseudo ranges, no atmospheric correction.
pub const MAX_SURVEY_BANCROFT_X_ERROR_M: f64 = 35.0;
pub const MAX_SURVEY_BANCROFT_Y_ERROR_M: f64 = 20.0;
pub const MAX_SURVEY_BANCROFT_Z_ERROR_M: f64 = 25.0;

/// [Config] used by the navigation tests on frozen data: default settings,
/// but the postfit smoother is disabled so every epoch evaluates the navigation
/// filter itself (the smoother freezes the published solution on its first value
/// with such sparse sampling).
pub fn navigation_test_config(method: crate::prelude::Method) -> crate::prelude::Config {
    let mut cfg = crate::prelude::Config::default().with_navigation_method(method);
    cfg.solver.postfit_denoising = 0.0;
    cfg
}

static INIT: Once = Once::new();

pub fn init_logger() {
    INIT.call_once(|| {
        env_logger::builder()
            .is_test(true)
            .filter_level(LevelFilter::Debug)
            .init();
    });
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
