//! PVT Solution
use crate::{
    navigation::{DilutionOfPrecision, State, sv::SVContribution},
    prelude::{Epoch, TimeScale},
};

#[cfg(feature = "serde")]
use serde::Serialize;

/// Describes the navigation technique used to obtain this [PVTSolution].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub enum PVTSolutionType {
    /// [PVTSolutionType::PPP] obtained using absolute navigation technique.
    /// This is true as long as this is not an RTK (differential) solution.
    /// In PPP solutions, the clock state, drift and TDOP are updated.
    PPP = 0,

    /// [PVTSolutionType::RTK] obtained using differential navigation technique.
    /// That means at minimum of one ground station was used to obtain
    /// this [PVTSolution]. The clock state, drift and TDoP are not resolved in
    /// RTK solutions. It is normal to obtain 0 for all these fields, for every single solution,
    /// for a session that only uses RTK.
    RTK = 1,
}

/// [PVTSolution] is solved by the navigation solver from a set of measurements,
/// using either PPP or RTK navigation technique. Depending on the technique being used,
/// the [PVTSolution] will differ. Mainly, RTK is not able to update and resolve the clock state.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize))]
pub struct PVTSolution {
    /// Type of solution
    pub solution_type: PVTSolutionType,

    /// Measurement [Epoch] that led to this solution.
    pub epoch: Epoch,

    /// Position solution, expressed in meters (ECEF).
    pub pos_m: (f64, f64, f64),

    /// Velocity solution, expressed in meters.s⁻¹ (ECEF).
    pub vel_m_s: (f64, f64, f64),

    /// Latitude, longitude and altitude above mean sea level,
    /// in degrees and meters.
    pub lat_long_alt_deg_deg_m: (f64, f64, f64),

    /// [TimeScale] of clock_offset [Duration] expression.
    pub timescale: TimeScale,

    /// Clock offset (in seconds).
    pub clock_offset_s: f64,

    /// Space Vehicles that helped form this solution
    /// and data associated to each individual SV
    pub sv: Vec<SVContribution>,

    /// Geometric Dilution of Precision
    pub gdop: f64,

    /// Vertical Dilution of Precision
    pub vdop: f64,

    /// Horizontal Diultion of Precision
    pub hdop: f64,

    /// Temporal Dilution of Precision
    pub tdop: f64,
}

impl PVTSolution {
    pub(crate) fn new(
        epoch: Epoch,
        uses_rtk: bool,
        state: &State,
        dop: &DilutionOfPrecision,
        contributions: &[SVContribution],
    ) -> Self {
        let pos_vel_ecef_m = state.to_position_velocity_ecef_m();
        let (clock_offset_s, _) = state.clock_profile_s();

        Self {
            epoch,
            solution_type: if uses_rtk {
                PVTSolutionType::RTK
            } else {
                PVTSolutionType::PPP
            },
            gdop: dop.gdop,
            tdop: dop.tdop,
            vdop: dop.vdop,
            hdop: dop.hdop,
            clock_offset_s,
            lat_long_alt_deg_deg_m: (
                state.lat_long_alt_deg_deg_km.0,
                state.lat_long_alt_deg_deg_km.1,
                state.lat_long_alt_deg_deg_km.2 * 1.0E3,
            ),
            sv: contributions.to_vec(),
            timescale: state.epoch.time_scale,
            pos_m: (pos_vel_ecef_m[0], pos_vel_ecef_m[1], pos_vel_ecef_m[2]),
            vel_m_s: (pos_vel_ecef_m[3], pos_vel_ecef_m[4], pos_vel_ecef_m[5]),
        }
    }
}
