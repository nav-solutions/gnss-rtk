//! Scenario description: everything the simulator needs to generate a run.
//!
//! A [Scenario] is a pure description (no state). The same scenario with the
//! same seed always produces bit identical measurements.
use std::str::FromStr;

use nalgebra::Vector3;

use crate::{
    prelude::{Carrier, Epoch, SV},
    tests::sim::{
        constellation::SimConstellation,
        geodesy::{ecef_to_enu_matrix, ecef_to_geodetic, geodetic_to_ecef},
    },
};

/// Receiver trajectory (truth)
#[derive(Debug, Clone, PartialEq)]
pub enum Trajectory {
    /// Static receiver, ECEF coordinates (m)
    Static { ecef_m: Vector3<f64> },

    /// Constant velocity, expressed in the local ENU frame of the starting point (m/s)
    Linear {
        start_ecef_m: Vector3<f64>,
        velocity_enu_m_s: Vector3<f64>,
    },

    /// Horizontal circle around a geodetic center (car on a ring road, pedestrian...)
    Circle {
        center_ecef_m: Vector3<f64>,
        radius_m: f64,
        period_s: f64,
    },
}

impl Trajectory {
    /// Static receiver at the given geodetic coordinates (deg, deg, m)
    pub fn static_geodetic(lat_deg: f64, lon_deg: f64, h_m: f64) -> Self {
        Self::Static {
            ecef_m: geodetic_to_ecef(lat_deg.to_radians(), lon_deg.to_radians(), h_m),
        }
    }

    /// Truth position (ECEF m), `t_s` seconds after scenario start.
    pub fn position_m(&self, t_s: f64) -> Vector3<f64> {
        match self {
            Self::Static { ecef_m } => *ecef_m,
            Self::Linear {
                start_ecef_m,
                velocity_enu_m_s,
            } => {
                let (lat, lon, _) = ecef_to_geodetic(start_ecef_m);
                let enu_to_ecef = ecef_to_enu_matrix(lat, lon).transpose();
                start_ecef_m + enu_to_ecef * (velocity_enu_m_s * t_s)
            },
            Self::Circle {
                center_ecef_m,
                radius_m,
                period_s,
            } => {
                let (lat, lon, _) = ecef_to_geodetic(center_ecef_m);
                let enu_to_ecef = ecef_to_enu_matrix(lat, lon).transpose();
                let phi = 2.0 * std::f64::consts::PI * t_s / period_s;
                let enu = Vector3::new(radius_m * phi.cos(), radius_m * phi.sin(), 0.0);
                center_ecef_m + enu_to_ecef * enu
            },
        }
    }

    /// Initial position (ECEF m)
    pub fn initial_position_m(&self) -> Vector3<f64> {
        self.position_m(0.0)
    }
}

/// Receiver clock truth: offset + drift.t + random walk
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClockModel {
    /// Initial offset (s)
    pub offset_s: f64,
    /// Drift (s/s)
    pub drift_s_s: f64,
    /// Random walk intensity (s / sqrt(s))
    pub random_walk_s_sqrt_s: f64,
}

impl ClockModel {
    /// Perfect clock
    pub fn perfect() -> Self {
        Self {
            offset_s: 0.0,
            drift_s_s: 0.0,
            random_walk_s_sqrt_s: 0.0,
        }
    }

    /// Typical steered TCXO: small offset, small drift, some random walk
    pub fn tcxo() -> Self {
        Self {
            offset_s: 150.0E-9,
            drift_s_s: 2.0E-10,
            random_walk_s_sqrt_s: 1.0E-10,
        }
    }
}

/// Measurement noise model: sigma(el) = sigma_zenith / sin(el) (when elevation dependent)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoiseModel {
    pub code_sigma_m: f64,
    pub phase_sigma_m: f64,
    pub elevation_dependent: bool,
}

impl NoiseModel {
    pub fn none() -> Self {
        Self {
            code_sigma_m: 0.0,
            phase_sigma_m: 0.0,
            elevation_dependent: false,
        }
    }

    /// Geodetic grade receiver
    pub fn geodetic() -> Self {
        Self {
            code_sigma_m: 0.3,
            phase_sigma_m: 0.003,
            elevation_dependent: true,
        }
    }

    /// Low cost receiver
    pub fn low_cost() -> Self {
        Self {
            code_sigma_m: 1.5,
            phase_sigma_m: 0.01,
            elevation_dependent: true,
        }
    }

    pub fn code_sigma(&self, elev_deg: f64) -> f64 {
        self.scale(elev_deg) * self.code_sigma_m
    }

    pub fn phase_sigma(&self, elev_deg: f64) -> f64 {
        self.scale(elev_deg) * self.phase_sigma_m
    }

    fn scale(&self, elev_deg: f64) -> f64 {
        if self.elevation_dependent {
            1.0 / elev_deg.to_radians().sin().max(0.1)
        } else {
            1.0
        }
    }
}

/// Ionosphere truth and model mismatch
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IonoModel {
    /// Vertical TEC (TECu)
    pub vtec_tecu: f64,
    /// Fraction of the truth that is NOT given to the solver through
    /// [EnvironmentalBias] (0: perfect model, 1: solver receives 0).
    pub model_error: f64,
}

impl IonoModel {
    pub fn none() -> Self {
        Self {
            vtec_tecu: 0.0,
            model_error: 0.0,
        }
    }

    /// Typical mid latitude day, perfectly modeled
    pub fn perfectly_modeled(vtec_tecu: f64) -> Self {
        Self {
            vtec_tecu,
            model_error: 0.0,
        }
    }

    /// Typical mid latitude day, the solver gets a model with 50 % error
    /// (Klobuchar like performance).
    pub fn klobuchar_like(vtec_tecu: f64) -> Self {
        Self {
            vtec_tecu,
            model_error: 0.5,
        }
    }
}

/// Troposphere: the solver receives the Niel model through [EnvironmentalBias],
/// the truth adds an unmodeled zenith wet delay.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TropoModel {
    /// Unmodeled zenith wet delay (m)
    pub residual_zwd_m: f64,
}

/// Cycle slip on both carriers of one satellite
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleSlip {
    pub sv: SV,
    /// Time of occurrence, seconds after scenario start
    pub at_s: f64,
    pub cycles_1: i64,
    pub cycles_2: i64,
}

/// Satellite signal outage (masking, tree, bridge...)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Outage {
    pub sv: SV,
    pub start_s: f64,
    pub end_s: f64,
}

/// Complete scenario description
#[derive(Debug, Clone)]
pub struct Scenario {
    pub name: String,
    /// Start of scenario (GPST)
    pub t0: Epoch,
    /// Duration (s)
    pub duration_s: f64,
    /// Sampling interval (s)
    pub interval_s: f64,
    /// Synthetic constellations
    pub constellations: Vec<SimConstellation>,
    /// Rover trajectory (truth)
    pub rover: Trajectory,
    /// Rover clock (truth)
    pub rover_clock: ClockModel,
    /// Optional static RTK base station (ECEF m)
    pub base: Option<Vector3<f64>>,
    /// Base clock (truth)
    pub base_clock: ClockModel,
    /// Noise model
    pub noise: NoiseModel,
    /// Ionosphere
    pub iono: IonoModel,
    /// Troposphere
    pub tropo: TropoModel,
    /// Simulator visibility mask (deg). The solver applies its own mask on top.
    pub elevation_mask_deg: f64,
    /// Cycle slips
    pub cycle_slips: Vec<CycleSlip>,
    /// Signal outages
    pub outages: Vec<Outage>,
    /// Satellite clock offset magnitude (s)
    pub sat_clock_offset_s: f64,
    /// Seed of all random processes
    pub seed: u64,
}

/// Default reference site (same area as the existing ESBC00DNK test data)
pub const REFERENCE_SITE_LAT_LON_H: (f64, f64, f64) = (55.4926, 8.4598, 60.0);

impl Scenario {
    /// Static, open sky, noise free, perfectly modeled, Galileo only.
    /// This is the "the math must be exact" baseline.
    pub fn ideal_static(name: &str) -> Self {
        let (lat, lon, h) = REFERENCE_SITE_LAT_LON_H;
        Self {
            name: name.to_string(),
            t0: Epoch::from_str("2024-03-14T00:00:00 GPST").unwrap(),
            duration_s: 3600.0,
            interval_s: 30.0,
            constellations: vec![SimConstellation::Galileo],
            rover: Trajectory::static_geodetic(lat, lon, h),
            rover_clock: ClockModel::perfect(),
            base: None,
            base_clock: ClockModel::perfect(),
            noise: NoiseModel::none(),
            iono: IonoModel::none(),
            tropo: TropoModel::default(),
            elevation_mask_deg: 5.0,
            cycle_slips: Vec::new(),
            outages: Vec::new(),
            sat_clock_offset_s: 5.0E-4,
            seed: 1,
        }
    }

    /// Static site with realistic geodetic receiver errors (noise, iono, clock).
    pub fn realistic_static(name: &str) -> Self {
        Self::ideal_static(name)
            .with_noise(NoiseModel::geodetic())
            .with_iono(IonoModel::klobuchar_like(20.0))
            .with_rover_clock(ClockModel::tcxo())
    }

    /// Adds a static RTK base at `east_m`, `north_m` from the rover initial position.
    pub fn with_base_offset_enu(mut self, east_m: f64, north_m: f64, up_m: f64) -> Self {
        let r0 = self.rover.initial_position_m();
        let (lat, lon, _) = ecef_to_geodetic(&r0);
        let enu_to_ecef = ecef_to_enu_matrix(lat, lon).transpose();
        self.base = Some(r0 + enu_to_ecef * Vector3::new(east_m, north_m, up_m));
        self
    }

    pub fn with_duration(mut self, duration_s: f64, interval_s: f64) -> Self {
        self.duration_s = duration_s;
        self.interval_s = interval_s;
        self
    }

    pub fn with_constellations(mut self, constellations: &[SimConstellation]) -> Self {
        self.constellations = constellations.to_vec();
        self
    }

    pub fn with_trajectory(mut self, trajectory: Trajectory) -> Self {
        self.rover = trajectory;
        self
    }

    pub fn with_rover_clock(mut self, clock: ClockModel) -> Self {
        self.rover_clock = clock;
        self
    }

    pub fn with_noise(mut self, noise: NoiseModel) -> Self {
        self.noise = noise;
        self
    }

    pub fn with_iono(mut self, iono: IonoModel) -> Self {
        self.iono = iono;
        self
    }

    pub fn with_tropo(mut self, tropo: TropoModel) -> Self {
        self.tropo = tropo;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    pub fn with_cycle_slip(mut self, slip: CycleSlip) -> Self {
        self.cycle_slips.push(slip);
        self
    }

    pub fn with_outage(mut self, outage: Outage) -> Self {
        self.outages.push(outage);
        self
    }

    /// Number of epochs
    pub fn num_epochs(&self) -> usize {
        (self.duration_s / self.interval_s).floor() as usize + 1
    }

    /// Epoch tag (receiver time) of the k-th sample
    pub fn epoch(&self, k: usize) -> Epoch {
        self.t0 + hifitime::Duration::from_seconds(k as f64 * self.interval_s)
    }

    /// Signals emitted by each constellation (primary, secondary)
    pub fn carriers(constellation: SimConstellation) -> (Carrier, Carrier) {
        match constellation {
            SimConstellation::Galileo => (Carrier::L1, Carrier::E5b),
            SimConstellation::GPS => (Carrier::L1, Carrier::L2),
        }
    }
}
