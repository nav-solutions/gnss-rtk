use crate::prelude::{Duration, IonosphereBias, SV, Signal};

#[cfg(feature = "serde")]
use serde::Serialize;

/// SV Navigation information
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub struct SVContribution {
    /// [SV] identity
    pub sv: SV,

    /// [Signal] being used
    pub signal: Signal,

    /// Orbital state in kilometers ECEF
    pub sv_pos_km: (f64, f64, f64),

    /// Orbital velocity in km/s ECEF
    pub sv_vel_km_s: (f64, f64, f64),

    /// Elevation angle from RX position
    pub elevation_deg: f64,

    /// Azimuth angle from RX position
    pub azimuth_deg: f64,

    /// Relativistic path range
    pub relativistic_path_range_m: f64,

    /// Troposphere bias as meters of propagation delay.
    pub tropo_bias: Option<f64>,

    /// Ionosphere bias as meters of propagation delay.
    pub iono_bias: Option<IonosphereBias>,

    /// Offset to selected [TimeScale], expressed as [Duration] within that [TimeScale].
    pub clock_correction: Option<Duration>,
}
