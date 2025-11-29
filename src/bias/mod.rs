/// Environmental biases
pub mod environment;

/// Spaceborn (onboard) biases
pub mod spaceborn;

use crate::prelude::{Epoch, SV, Vector3};

/// [BiasRuntime] contains everything to compute your [Bias] estimate.
#[derive(Default, Copy, Clone)]
pub struct BiasRuntime {
    /// [SV] identity
    pub sv: SV,

    /// [Epoch] of measurement.
    pub epoch: Epoch,

    /// Selected [SV] position, in meters ECEF.
    pub sv_position_ecef_m: (f64, f64, f64),

    /// Selected [SV] (elevation, azimuth) attitude, both in degrees.
    pub sv_elevation_azimuth_deg_deg: (f64, f64),

    /// Estimated receiver state, in meters ECEF.
    pub rx_position_m: Vector3<f64>,

    /// Latest estimate of the receiver (latitude, longitude, altitude above mean sea level),
    /// coordinates, in degrees, degrees and kilometers.
    pub rcvr_lat_long_alt_deg_deg_km: (f64, f64, f64),

    /// Signal frequency used in the navigation process, in Hertz.
    pub frequency_hz: f64,
}
