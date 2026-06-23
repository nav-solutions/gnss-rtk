use crate::cfg::Error;
use crate::prelude::{Duration, SPEED_OF_LIGHT_M_S};

use std::f64::consts::PI;

use nalgebra::DMatrix;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Default acceleration PSD
const fn default_accel_psd() -> f64 {
    UserProfile::Pedestrian.psd()
}

/// Default clock PSD
const fn default_clock_psd() -> f64 {
    ClockProfile::Quartz.bias_psd()
}

/// Default clock PSD
const fn default_clock_drift_psd() -> f64 {
    ClockProfile::Quartz.drift_psd()
}

/// [UserProfile] can be used to generate a set of [UserParameters] easily.
#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum UserProfile {
    /// [UserProfile::Static] applies to user applications where
    /// the receiver antenna is held static at all times.
    Static,

    /// [UserProfile::Pedestrian]: < 10 km/h very low velocity.
    /// This is our default mode.
    #[cfg_attr(feature = "serde", serde(alias = "pedestrian", alias = "Pedestrian"))]
    #[default]
    Pedestrian,

    /// [UserProfile::Car]: < 100 km/h low velocity
    #[cfg_attr(feature = "serde", serde(alias = "car", alias = "Car"))]
    Car,

    /// [UserProfile::Airplane]: < 1000 km/h high velocity
    #[cfg_attr(feature = "serde", serde(alias = "airplane", alias = "airplane"))]
    Airplane,

    /// [UserProfile::Rocket]: > 1000 km/h ultra high velocity
    #[cfg_attr(feature = "serde", serde(alias = "rocket", alias = "rocket"))]
    Rocket,
}

impl UserProfile {
    pub(crate) const fn psd(&self) -> f64 {
        match self {
            Self::Static => 1.0,
            Self::Pedestrian => 2.0,
            Self::Car => 10.0,
            Self::Airplane => 100.0,
            Self::Rocket => 1_000.0,
        }
    }
}

impl std::str::FromStr for UserProfile {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_lowercase();
        let trimmed = s.trim();
        match trimmed {
            "static" => Ok(Self::Static),
            "pedestrian" => Ok(Self::Pedestrian),
            "car" => Ok(Self::Car),
            "airplane" => Ok(Self::Airplane),
            "rocket" => Ok(Self::Rocket),
            _ => Err(Error::InvalidUserProfile),
        }
    }
}

impl std::fmt::Display for UserProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Static => write!(f, "static"),
            Self::Pedestrian => write!(f, "pedestrian"),
            Self::Car => write!(f, "car"),
            Self::Airplane => write!(f, "airplane"),
            Self::Rocket => write!(f, "rocket"),
        }
    }
}

#[derive(Default, Debug, Clone, Copy, PartialEq)]
pub enum ClockProfile {
    /// [ClockProfile::Quartz] low quality clock,
    /// poorer than GNSS constellation.
    #[default]
    Quartz,

    /// [ClockProfile::Oscillator] medium quality clock,
    /// poorer than GNSS constellation.
    Oscillator,

    /// [ClockProfile::Atomic] high quality clock, potentially
    /// at the same level of a GNSS constellation.
    Atomic,

    /// [ClockProfile::HydrogenMaser] ultra high quality clock,
    /// better than a GNSS constellation.
    HydrogenMaser,
}

impl std::str::FromStr for ClockProfile {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_lowercase();
        let trimmed = s.trim();
        match trimmed {
            "quartz" => Ok(Self::Quartz),
            "maser" | "h-maser" => Ok(Self::HydrogenMaser),
            "oscillator" | "ocxo" => Ok(Self::Oscillator),
            "atomic" | "rb" | "rubidium" => Ok(Self::Atomic),
            _ => Err(Error::InvalidClockProfile),
        }
    }
}

impl std::fmt::Display for ClockProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Quartz => write!(f, "quartz"),
            Self::Oscillator => write!(f, "oscillator"),
            Self::Atomic => write!(f, "atomic"),
            Self::HydrogenMaser => write!(f, "hydrogen maser"),
        }
    }
}

impl ClockProfile {
    /// Returns bias PSD for this [ClockProfile].
    pub(crate) const fn bias_psd(&self) -> f64 {
        match self {
            Self::Quartz => 0.5 * 2.0 * 2.0E-19,
            Self::Oscillator => 0.5 * 8.0 * 2.0E-20,
            Self::Atomic => 0.5 * 2.0 * 2.0E-20,
            Self::HydrogenMaser => 0.5 * 2.0E-20,
        }
    }

    /// Returns drift PSD for this [ClockProfile].
    pub(crate) const fn drift_psd(&self) -> f64 {
        match self {
            Self::Quartz => 2.0 * PI * PI * 2.0E-20,
            Self::Oscillator => 2.0 * PI * PI * 4.0 * 2.0E-23,
            Self::Atomic => 2.0 * PI * PI * 4.0 * 2.0E-29,
            Self::HydrogenMaser => 2.0 * PI * PI * 4.0 * 2.0E-29,
        }
    }
}

/// [UserParameters] are applications dependent.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct UserParameters {
    /// Acceleration Power Spectral Density (PSD) in m².s⁻.¹
    /// Directly tied to your application profile and attitude.
    /// For cars we recommend 2.0. For pedestrian 0.5,
    /// and obviously 0.0 for static applications.
    #[cfg_attr(feature = "serde", serde(default = "default_accel_psd"))]
    pub accel_psd: f64,

    /// Clock bias PSD
    #[cfg_attr(feature = "serde", serde(default = "default_clock_psd"))]
    pub clock_psd: f64,

    /// Clock drift PSD
    #[cfg_attr(feature = "serde", serde(default = "default_clock_drift_psd"))]
    pub clock_drift_psd: f64,
}

impl std::fmt::Display for UserParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "a={}m².s⁻¹, offset={}s, drift={}s.s⁻¹",
            self.accel_psd, self.clock_psd, self.clock_drift_psd
        )
    }
}

impl Default for UserParameters {
    fn default() -> Self {
        Self {
            accel_psd: default_accel_psd(),
            clock_psd: default_clock_psd(),
            clock_drift_psd: default_clock_drift_psd(),
        }
    }
}

impl UserParameters {
    /// Creates a [UserParameters] set from [UserProfile] and [ClockProfile]
    pub fn new(user_profile: UserProfile, clock_profile: ClockProfile) -> Self {
        Self {
            accel_psd: user_profile.psd(),
            clock_psd: clock_profile.bias_psd(),
            clock_drift_psd: clock_profile.drift_psd(),
        }
    }

    /// Parametrization of the coveriance [DMatrix]
    pub(crate) fn q_matrix(&self, q_mat: &mut DMatrix<f64>, dt: Duration, ndf: usize) {
        assert!(ndf > 2, "Q cov: minimal dimension");

        let dt_s = dt.to_seconds();
        let dt_s3 = dt_s.powi(3);

        for i in 0..=2 {
            q_mat[(i, i)] = self.accel_psd;

            if self.accel_psd != 1.0 {
                q_mat[(i, i)] *= dt_s3 / 3.0;
            }
        }

        if ndf > 3 {
            // if dt == Duration::ZERO {
            //    q_mat[(3, 3)] = (100.0e-3 * SPEED_OF_LIGHT_M_S).powi(2);
            //} else {
            q_mat[(3, 3)] = SPEED_OF_LIGHT_M_S.powi(2)
                * (self.clock_psd * dt_s + self.clock_drift_psd * dt_s3 / 3.0);
            //}
        }
    }
}
