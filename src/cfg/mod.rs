use thiserror::Error;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

mod method;
mod modeling;
mod solver;

pub use crate::{
    carrier::Signal,
    cfg::solver::{SolverOpts, Weighting},
    cfg::{method::Method, modeling::Modeling},
    prelude::TimeScale,
};

/// Configuration Error
#[derive(Debug, Error)]
pub enum Error {
    #[error("invalid troposphere model")]
    InvalidTroposphereModel,

    #[error("invalid user profile")]
    InvalidUserProfile,

    #[error("invalid clock profile")]
    InvalidClockProfile,
}

const fn default_timescale() -> TimeScale {
    TimeScale::GPST
}

const fn max_tropo_bias() -> f64 {
    30.0
}

const fn max_iono_bias() -> f64 {
    10.0
}

const fn min_sv_elev() -> Option<f64> {
    Some(12.5)
}

const fn default_code_smoothing() -> usize {
    0
}

const fn default_eclipse_rate_percent() -> f64 {
    10.0
}

#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
/// System Internal Delay as defined by BIPM in
/// "GPS Receivers Accurate Time Comparison" : the (frequency dependent)
/// time delay introduced by the combination of:
///  + the RF cable (up to several nanoseconds)
///  + the distance between the antenna baseline and its APC:
///    a couple picoseconds, and is frequency dependent
///  + the GNSS receiver inner delay (hardware and frequency dependent)
pub struct InternalDelay {
    /// Delay [s]
    pub delay: f64,

    /// Carrier frequency [Hz]
    pub frequency: f64,
}

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Config {
    /// Time scale in which we express the PVT solutions,
    /// [TimeScale::GPST] is the default value.
    #[cfg_attr(feature = "serde", serde(default = "default_timescale"))]
    pub timescale: TimeScale,

    /// Navigation [Method] (technique) to be used.
    #[cfg_attr(feature = "serde", serde(default))]
    pub method: Method,

    /// Select a prefered signal.
    /// When defined, this signal will strictly be used in the navigation process.
    /// When undefined, the algorithm will prefer the best SNR available, and the
    /// signal frequency being used might change.
    /// When [Method] is [Method::SPP] this should be a single frequency.
    /// When [Method] is not [Method::SPP] and [Modeling] enables Ionospheric
    /// bias compensation,
    #[cfg_attr(feature = "serde", serde(default))]
    pub prefered_signal: Option<Signal>,

    /// Fixed altitude: reduces the need of 4 to 3 SV to obtain 3D solutions.
    #[cfg_attr(feature = "serde", serde(default))]
    pub fixed_altitude: Option<f64>,

    /// Pseudo Range code smoothing (window length).
    /// Use phase observatoins to smooth and reduce error in the pseudo range code.
    /// This has no effect if phase observations are missing.
    /// Set to 0 to disable this feature completely.
    /// When parametrizing, think in terms of window duration versus Ionospheric activity.
    #[cfg_attr(feature = "serde", serde(default = "default_code_smoothing"))]
    pub code_smoothing: usize,

    /// Internal delays to compensate for (total summation, in [s]).
    /// Compensation is only effective if [Modeling.cable_delay]
    /// is also turned on.
    #[cfg_attr(feature = "serde", serde(default))]
    pub int_delay: Vec<InternalDelay>,

    /// Antenna Reference Point (ARP) expressed as ENU offset [m]
    #[cfg_attr(feature = "serde", serde(default))]
    pub arp_enu: Option<(f64, f64, f64)>,

    /// Solver customization
    #[cfg_attr(feature = "serde", serde(default))]
    pub solver: SolverOpts,

    /// Time Reference Delay. According to BIPM ""GPS Receivers Accurate Time Comparison""
    /// this is the time delay between the receiver external reference clock
    /// and the internal sampling clock. This is typically needed in
    /// precision timing applications.
    #[cfg_attr(feature = "serde", serde(default))]
    pub externalref_delay_s: Option<f64>,

    /// Maximal Earth / Sun occultation tolerated for each satellite.
    /// For example, 20.0% means that we require satellites to be 80% illmuinated.
    /// 10.0% is our default value.
    #[cfg_attr(
        feature = "serde",
        serde(alias = "max_eclipse_rate", default = "default_eclipse_rate_percent")
    )]
    pub max_eclipse_rate_percent: f64,

    /// Minimal SV elevation angle for an SV to contribute to the solution.
    /// Use this as a simple quality criteria.
    #[cfg_attr(feature = "serde", serde(default = "min_sv_elev"))]
    pub min_sv_elev: Option<f64>,

    /// Minimal SV Azimuth angle for an SV to contribute to the solution.
    /// SV below that angle will not be considered.
    /// Use this is in special navigation scenarios.
    #[cfg_attr(feature = "serde", serde(default))]
    pub min_sv_azim: Option<f64>,

    /// Maximal SV Azimuth angle for an SV to contribute to the solution.
    /// SV below that angle will not be considered.
    /// Use this is in special navigation scenarios.
    #[cfg_attr(feature = "serde", serde(default))]
    pub max_sv_azim: Option<f64>,

    /// Minimal SNR for an SV to contribute to the solution.
    #[cfg_attr(feature = "serde", serde(default))]
    pub min_snr: Option<f64>,

    /// Maximal tropo bias that we tolerate (in [m]).
    /// Has no effect if modeling.tropo_delay is disabled.
    #[cfg_attr(feature = "serde", serde(default = "max_tropo_bias"))]
    pub max_tropo_bias: f64,

    /// Maximal iono bias that we tolerate (in [m]).
    /// Has no effect if modeling.iono_delay is disabled.
    #[cfg_attr(feature = "serde", serde(default = "max_iono_bias"))]
    pub max_iono_bias: f64,

    /// Atmospherical and Physical [Modeling] used to improve the accuracy of solution.
    #[cfg_attr(feature = "serde", serde(default))]
    pub modeling: Modeling,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            timescale: default_timescale(),
            method: Method::default(),
            solver: SolverOpts::default(),
            int_delay: Default::default(),
            modeling: Modeling::default(),
            fixed_altitude: None,
            prefered_signal: None,
            arp_enu: None,
            min_snr: None, // TODO
            min_sv_azim: None,
            max_sv_azim: None,
            externalref_delay_s: None,
            min_sv_elev: min_sv_elev(),
            max_iono_bias: max_iono_bias(),
            max_tropo_bias: max_tropo_bias(),
            code_smoothing: default_code_smoothing(),
            max_eclipse_rate_percent: default_eclipse_rate_percent(),
        }
    }
}

impl Config {
    /// Returns new [Config] with desired navigation [Method]
    pub fn with_navigation_method(&self, method: Method) -> Self {
        let mut s = self.clone();
        s.method = method;
        s
    }

    /// Copies and returns [Config] with desired [Modeling], that you can
    /// tune to improve the accuracy of your solution, or for learning purposes.
    pub fn with_modeling(&self, modeling: Modeling) -> Self {
        let mut s = self.clone();
        s.modeling = modeling;
        s
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod test {
    use super::*;
    use std::io::Write;

    #[test]
    fn generate_default_preset() {
        let cfg = Config::default();
        let string = serde_json::to_string_pretty(&cfg).unwrap();
        let mut fd = std::fs::File::create("default.json").unwrap();
        write!(fd, "{}", string).unwrap();
    }
}
