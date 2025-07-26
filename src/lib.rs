#![doc(
    html_logo_url = "https://raw.githubusercontent.com/nav-solutions/.github/master/logos/logo2.jpg"
)]
#![doc = include_str!("../README.md")]
#![cfg_attr(docsrs, feature(doc_cfg))]

extern crate gnss_rs as gnss;

pub mod error;

// mod ambiguity;
// mod averager;
mod bancroft;
mod bias;
mod candidate;
mod carrier;
mod cfg;
mod ephemeris;
mod navigation;
mod orbit;
mod pool;
mod rtk;
// mod smoothing;
mod time;
mod user;
// mod tides;

pub(crate) mod constants;
pub(crate) mod solver;

#[cfg(test)]
mod tests;

// prelude
pub mod prelude {
    pub use crate::{
        bias::{
            environment::{
                EnvironmentalBias, IonosphereBias, IonosphereModel, KbModel, TroposphereModel,
            },
            spaceborn::{SatelliteClockCorrection, SpacebornBias},
            BiasRuntime,
        },
        candidate::{Candidate, Observation},
        carrier::{Carrier, Signal},
        cfg::{Config, Method},
        constants::SPEED_OF_LIGHT_M_S,
        ephemeris::{Ephemeris, EphemerisSource},
        error::Error,
        navigation::solutions::{PVTSolution, PVTSolutionType},
        orbit::OrbitSource,
        rtk::RTKBase,
        solver::Solver,
        time::AbsoluteTime,
        user::{ClockProfile, UserParameters, UserProfile},
    };

    // std types
    pub use std::rc::Rc;

    // gnss types
    pub use gnss::prelude::{Constellation, SV};

    // anise types
    pub use anise::{
        constants::frames::{EARTH_ITRF93, EARTH_J2000, IAU_EARTH_FRAME, SUN_J2000},
        naif::SPK,
        prelude::{Aberration, Almanac, Frame, Orbit},
    };

    // hifitime types
    pub use hifitime::{Duration, Epoch, TimeScale};

    // nalgebra
    pub use nalgebra::{Vector3, Vector4, Vector6};
}
