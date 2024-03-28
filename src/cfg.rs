use thiserror::Error;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::prelude::TimeScale;
use nalgebra::DMatrix;

/// Configuration Error
#[derive(Debug, Error)]
pub enum Error {
    #[error("unknown tropo model")]
    UnknownTropoModel(String),
}

/// Solving method
#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum Method {
    /// Single Point Positioning (SPP).
    /// Uses pseudorange observations from a single carrier signal.
    /// There is not point providing carrier phase obserations with this method.
    /// SPP can work in degraged environments where a unique signal is sampled.
    /// Combine this to advanced configurations and you may obtain metric precision.
    #[default]
    SPP,
    /// Precise Point Positioning (PPP)
    PPP,
}

impl std::fmt::Display for Method {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::SPP => write!(fmt, "SPP"),
            Self::PPP => write!(fmt, "PPP"),
        }
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
/// Positioning method
pub enum Positioning {
    /// Receiver is static
    #[default]
    Static,
    /// Receiver is moving
    Kinematic,
}

/// Filter to use in the solving process
#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum Filter {
    /// No filter. The "filter" part of the
    /// [Config] struct is disregarded
    None,
    /// LSQ Filter. Heavy computation, converges much slower than a Kalman filter.
    #[default]
    LSQ,
    /// Kalman Filter. Heavy+ computation, converges much faster than LSQ.
    KF,
}

#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub struct ElevationMappingFunction {
    /// a + b * e-elev/c
    pub a: f64,
    /// a + b * e-elev/c
    pub b: f64,
    /// a + b * e-elev/c
    pub c: f64,
}

impl ElevationMappingFunction {
    pub(crate) fn eval(&self, elev_sv: f64) -> f64 {
        self.a + self.b * (elev_sv / self.c).exp()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum WeightMatrix {
    /// a + b e-elev/c
    MappingFunction(ElevationMappingFunction),
    /// Advanced measurement noise covariance matrix
    Covar,
}

fn default_timescale() -> TimeScale {
    TimeScale::GPST
}

fn default_interp() -> usize {
    11
}

fn default_max_sv() -> usize {
    10
}

fn default_smoothing() -> bool {
    false
}

fn default_sv_clock() -> bool {
    true
}

fn default_sv_tgd() -> bool {
    true
}

fn default_iono() -> bool {
    true
}

fn default_tropo() -> bool {
    true
}

fn default_earth_rot() -> bool {
    true
}

fn default_relativistic_clock_bias() -> bool {
    true
}

fn default_relativistic_path_range() -> bool {
    false
}

fn default_sv_apc() -> bool {
    false
}

fn default_weight_matrix() -> Option<WeightMatrix> {
    None
    //Some(WeightMatrix::MappingFunction(
    //    ElevationMappingFunction {
    //        a: 5.0,
    //        b: 0.0,
    //        c: 10.0,
    //    },
    //))
}

fn default_filter_opts() -> Option<FilterOpts> {
    Some(FilterOpts {
        weight_matrix: default_weight_matrix(),
    })
}

fn default_gdop_threshold() -> Option<f64> {
    None
}

fn default_tdop_threshold() -> Option<f64> {
    None
}

#[derive(Default, Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
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

#[derive(Default, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub struct SolverOpts {
    /// GDOP threshold to invalidate ongoing GDOP
    #[cfg_attr(feature = "serde", serde(default = "default_gdop_threshold"))]
    pub gdop_threshold: Option<f64>,
    /// TDOP threshold to invalidate ongoing TDOP
    #[cfg_attr(feature = "serde", serde(default = "default_tdop_threshold"))]
    pub tdop_threshold: Option<f64>,
    /// Filter to use
    #[cfg_attr(feature = "serde", serde(default))]
    pub filter: Filter,
    /// Filter options
    #[cfg_attr(feature = "serde", serde(default = "default_filter_opts"))]
    pub filter_opts: Option<FilterOpts>,
}

#[derive(Default, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub struct FilterOpts {
    /// Weight Matrix
    #[cfg_attr(feature = "serde", serde(default = "default_weight_matrix"))]
    pub weight_matrix: Option<WeightMatrix>,
}

impl SolverOpts {
    /*
     * form the weight matrix to be used in the solving process
     */
    pub(crate) fn weight_matrix(&self, _nb_rows: usize, sv_elev: Vec<f64>) -> DMatrix<f64> {
        let mut mat = DMatrix::identity(sv_elev.len(), sv_elev.len());
        if let Some(opts) = &self.filter_opts {
            match &opts.weight_matrix {
                Some(WeightMatrix::Covar) => panic!("not implemented yet"),
                Some(WeightMatrix::MappingFunction(mapf)) => {
                    for i in 0..sv_elev.len() - 1 {
                        let sigma = mapf.a + mapf.b * ((-sv_elev[i]) / mapf.c).exp();
                        mat[(i, i)] = 1.0 / sigma.powi(2);
                    }
                },
                None => {},
            }
        }
        mat
    }
}

/// Atmospherical, Physical and Environmental modeling
#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Modeling {
    #[cfg_attr(feature = "serde", serde(default))]
    pub sv_clock_bias: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub sv_total_group_delay: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub sv_apc: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub relativistic_clock_bias: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub relativistic_path_range: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub tropo_delay: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub iono_delay: bool,
    #[cfg_attr(feature = "serde", serde(default))]
    pub earth_rotation: bool,
}

impl Default for Modeling {
    fn default() -> Self {
        Self {
            sv_clock_bias: default_sv_clock(),
            sv_apc: default_sv_apc(),
            iono_delay: default_iono(),
            tropo_delay: default_tropo(),
            sv_total_group_delay: default_sv_tgd(),
            earth_rotation: default_earth_rot(),
            relativistic_clock_bias: default_relativistic_clock_bias(),
            relativistic_path_range: default_relativistic_path_range(),
        }
    }
}

impl Modeling {
    pub fn preset(_method: Method, _filter: Filter) -> Self {
        Self::default()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub struct Config {
    /// Time scale
    #[cfg_attr(feature = "serde", serde(default = "default_timescale"))]
    pub timescale: TimeScale,
    /// Method to use
    #[cfg_attr(feature = "serde", serde(default))]
    pub method: Method,
    /// (Position) interpolation filter order.
    /// A minimal order must be respected for correct results.
    /// -  7 is the minimal value for metric resolution
    /// - 11 is the standard when aiming at submetric resolution
    #[cfg_attr(feature = "serde", serde(default = "default_interp"))]
    pub interp_order: usize,
    /// Fixed altitude : reduces the need of 4 to 3 SV to resolve a solution
    #[cfg_attr(feature = "serde", serde(default))]
    pub fixed_altitude: Option<f64>,
    /// PR code smoothing filter before moving forward
    #[cfg_attr(feature = "serde", serde(default = "default_smoothing"))]
    pub code_smoothing: bool,
    /// Internal delays
    #[cfg_attr(feature = "serde", serde(default))]
    pub int_delay: Vec<InternalDelay>,
    /// Antenna Reference Point (ARP) as ENU offset [m]
    #[cfg_attr(feature = "serde", serde(default))]
    pub arp_enu: Option<(f64, f64, f64)>,
    /// Solver customization
    #[cfg_attr(feature = "serde", serde(default))]
    pub solver: SolverOpts,
    /// Time Reference Delay, as defined by BIPM in
    /// "GPS Receivers Accurate Time Comparison" : the time delay
    /// between the GNSS receiver external reference clock and the sampling clock
    /// (once again can be persued as a cable delay). This one is typically
    /// only required in ultra high precision timing applications
    #[cfg_attr(feature = "serde", serde(default))]
    pub externalref_delay: Option<f64>,
    /// Minimal percentage ]0; 1[ of Sun light to be received by an SV
    /// for not to be considered in Eclipse.
    /// A value closer to 0 means we tolerate fast Eclipse exit.
    /// A value closer to 1 is a stringent criteria: eclipse must be totally exited.
    #[cfg_attr(feature = "serde", serde(default))]
    pub min_sv_sunlight_rate: Option<f64>,
    /// Minimal elevation angle. SV below that angle will not be considered.
    #[cfg_attr(feature = "serde", serde(default))]
    pub min_sv_elev: Option<f64>,
    /// Minimal SNR for an SV to be considered.
    #[cfg_attr(feature = "serde", serde(default))]
    pub min_snr: Option<f64>,
    /// modeling
    #[cfg_attr(feature = "serde", serde(default))]
    pub modeling: Modeling,
    /// Max. number of vehicules to consider.
    /// The more the merrier, but it also means heavier computations
    #[cfg_attr(feature = "serde", serde(default = "default_max_sv"))]
    pub max_sv: usize,
}

impl Config {
    pub fn preset(method: Method) -> Self {
        match method {
            Method::SPP => Self {
                method,
                timescale: default_timescale(),
                fixed_altitude: None,
                interp_order: default_interp(),
                code_smoothing: default_smoothing(),
                min_sv_sunlight_rate: None,
                min_sv_elev: Some(15.0),
                min_snr: Some(30.0),
                modeling: Modeling::default(),
                max_sv: default_max_sv(),
                int_delay: Default::default(),
                externalref_delay: Default::default(),
                arp_enu: None,
                solver: SolverOpts {
                    gdop_threshold: default_gdop_threshold(),
                    tdop_threshold: default_tdop_threshold(),
                    filter: Filter::LSQ,
                    filter_opts: default_filter_opts(),
                },
            },
            Method::PPP => panic!("not available yet"),
        }
    }
}
