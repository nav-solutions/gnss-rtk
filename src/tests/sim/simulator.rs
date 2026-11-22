//! Measurement simulator: forward GNSS observation model with known truth.
//!
//! The simulator implements every input trait of the [Solver]
//! ([OrbitSource], [SpacebornBias], [EnvironmentalBias], [RTKBase]), so the
//! solver is exercised exactly like in a real application, through its
//! public API only.
//!
//! Forward model, for station `r`, satellite `s`, carrier `f`, at receiver time tag `T`:
//!
//! ```text
//! t_rx  = T - dt_r                         (true reception time)
//! tau   = |R3(we.tau) x_s(t_rx - tau) - x_r(t_rx)| / c   (iterated, Sagnac included)
//! P_f   = c.tau + S + c.(dt_r - dt_s) + T_trop + I_f + e_P
//! L_f   = c.tau + S + c.(dt_r - dt_s) + T_trop - I_f + lambda_f.(N_f + slips) + e_L
//! ```
//!
//! with S the relativistic (Shapiro) path delay.
//!
//! ```text
//! ```
use std::{cell::RefCell, collections::HashMap};

use nalgebra::Vector3;

use crate::{
    constants::{EARTH_ANGULAR_VEL_RAD, SPEED_OF_LIGHT_M_S},
    prelude::{
        BiasRuntime, Candidate, Carrier, Duration, EnvironmentalBias, Epoch, Frame, Observation,
        Orbit, OrbitSource, RTKBase, Rc, SV, SatelliteClockCorrection, SpacebornBias,
        TroposphereModel,
    },
    tests::sim::{
        constellation::CircularOrbit,
        geodesy::{ecef_to_geodetic, elevation_azimuth_deg},
        rng::SimRng,
        scenario::{ClockModel, Scenario},
    },
};

/// Earth gravitational parameter (m³.s⁻²), truth side
const GM_M3_S2: f64 = 3.986_004_418E14;

/// Ionosphere thin shell height (m)
const IONO_SHELL_HEIGHT_M: f64 = 350_000.0;
const EARTH_MEAN_RADIUS_M: f64 = 6_371_000.0;

/// Measurement station
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Station {
    Rover,
    Base,
}

impl Station {
    fn salt(&self) -> u64 {
        match self {
            Self::Rover => 0x1111,
            Self::Base => 0x2222,
        }
    }
}

/// Truth associated to one simulated epoch
#[derive(Debug, Clone)]
pub struct Truth {
    /// Receiver position at true reception time (ECEF m)
    pub position_m: Vector3<f64>,
    /// Receiver clock offset (s)
    pub clock_offset_s: f64,
}

/// One simulated epoch
#[derive(Debug, Clone)]
pub struct SimEpoch {
    /// Sample index
    pub k: usize,
    /// Receiver time tag
    pub epoch: Epoch,
    /// Seconds since scenario start
    pub t_s: f64,
    /// Rover measurements
    pub candidates: Vec<Candidate>,
    /// Rover truth
    pub truth: Truth,
}

/// Per satellite / per signal truth, exposed for white box tests.
#[derive(Debug, Clone, Copy)]
pub struct SignalTruth {
    pub sv: SV,
    pub geometric_range_m: f64,
    /// Satellite position at transmission time, rotated to the earth fixed
    /// frame at reception time (ECEF m)
    pub sat_position_m: Vector3<f64>,
    /// Relativistic (Shapiro) path delay (m)
    pub shapiro_m: f64,
    pub time_of_flight_s: f64,
    pub elevation_deg: f64,
    pub sat_clock_s: f64,
    pub tropo_m: f64,
    pub iono_l1_m: f64,
    pub ambiguity_1: i64,
    pub ambiguity_2: i64,
}

/// The simulator. Cheap to clone (shared state behind [Rc]).
#[derive(Clone)]
pub struct Simulator {
    inner: Rc<Inner>,
}

struct Inner {
    scenario: Scenario,
    frame: Frame,
    orbits: Vec<CircularOrbit>,
    constellation_of: HashMap<SV, (Carrier, Carrier)>,
    sat_clock: HashMap<SV, (f64, f64)>,
    rover_clock: Vec<f64>,
    base_clock: Vec<f64>,
    /// statistics: how many times the solver requested orbits
    orbit_requests: RefCell<usize>,
}

fn simulate_clock(model: &ClockModel, scenario: &Scenario, rng: &mut SimRng) -> Vec<f64> {
    let n = scenario.num_epochs();
    let mut out = Vec::with_capacity(n);
    let mut walk = 0.0;

    for k in 0..n {
        let t = k as f64 * scenario.interval_s;
        if k > 0 {
            walk += rng.normal(model.random_walk_s_sqrt_s * scenario.interval_s.sqrt());
        }
        out.push(model.offset_s + model.drift_s_s * t + walk);
    }

    out
}

impl Simulator {
    pub fn new(scenario: Scenario, frame: Frame) -> Self {
        let rng = SimRng::new(scenario.seed);

        let mut orbits = Vec::new();
        let mut constellation_of = HashMap::new();
        let mut sat_clock = HashMap::new();

        for constellation in scenario.constellations.iter() {
            let carriers = Scenario::carriers(*constellation);

            for orbit in constellation.orbits() {
                let mut sv_rng = rng.fork(0xC10C + orbit.sv.prn as u64 * 131 + carriers.1 as u64);

                let a0 =
                    sv_rng.uniform_in(-scenario.sat_clock_offset_s, scenario.sat_clock_offset_s);
                let a1 = sv_rng.uniform_in(-1.0E-12, 1.0E-12);

                sat_clock.insert(orbit.sv, (a0, a1));
                constellation_of.insert(orbit.sv, carriers);
                orbits.push(orbit);
            }
        }

        let rover_clock = simulate_clock(&scenario.rover_clock, &scenario, &mut rng.fork(0xAAAA));
        let base_clock = simulate_clock(&scenario.base_clock, &scenario, &mut rng.fork(0xBBBB));

        Self {
            inner: Rc::new(Inner {
                scenario,
                frame,
                orbits,
                constellation_of,
                sat_clock,
                rover_clock,
                base_clock,
                orbit_requests: RefCell::new(0),
            }),
        }
    }

    pub fn scenario(&self) -> &Scenario {
        &self.inner.scenario
    }

    pub fn frame(&self) -> Frame {
        self.inner.frame
    }

    /// Number of orbital states the solver requested so far
    pub fn orbit_requests(&self) -> usize {
        *self.inner.orbit_requests.borrow()
    }

    fn seconds_since_start(&self, t: Epoch) -> f64 {
        (t - self.inner.scenario.t0).to_seconds()
    }

    fn index_of(&self, t: Epoch) -> usize {
        let k = (self.seconds_since_start(t) / self.inner.scenario.interval_s).round();
        (k.max(0.0) as usize).min(self.inner.scenario.num_epochs() - 1)
    }

    fn orbit(&self, sv: SV) -> Option<&CircularOrbit> {
        self.inner.orbits.iter().find(|orb| orb.sv == sv)
    }

    /// Satellite clock offset (s) at receiver time tag `t`, as provided to the solver
    /// (full precision: precise clock products resolve 1E-12 s, far below the
    /// 1 ns resolution of [Duration]).
    pub fn sat_clock_s(&self, sv: SV, t: Epoch) -> f64 {
        let (a0, a1) = self.inner.sat_clock.get(&sv).copied().unwrap_or_default();
        let dt = self.seconds_since_start(t);
        a0 + a1 * dt
    }

    /// Truth position of a station, `t_s` seconds after start
    pub fn station_position_m(&self, station: Station, t_s: f64) -> Vector3<f64> {
        match station {
            Station::Rover => self.inner.scenario.rover.position_m(t_s),
            Station::Base => self
                .inner
                .scenario
                .base
                .expect("scenario does not define a base station"),
        }
    }

    /// Truth clock of a station at sample `k`
    pub fn station_clock_s(&self, station: Station, k: usize) -> f64 {
        match station {
            Station::Rover => self.inner.rover_clock[k],
            Station::Base => self.inner.base_clock[k],
        }
    }

    /// Ionospheric delay (m) for given elevation and frequency
    pub fn iono_delay_m(&self, elev_deg: f64, frequency_hz: f64) -> f64 {
        let vtec = self.inner.scenario.iono.vtec_tecu;
        if vtec == 0.0 {
            return 0.0;
        }
        let ratio = EARTH_MEAN_RADIUS_M / (EARTH_MEAN_RADIUS_M + IONO_SHELL_HEIGHT_M);
        let obliquity = 1.0 / (1.0 - (ratio * elev_deg.to_radians().cos()).powi(2)).sqrt();
        40.3E16 * vtec * obliquity / frequency_hz.powi(2)
    }

    /// Tropospheric delay truth (m): Niel model + unmodeled wet residual
    pub fn tropo_delay_m(
        &self,
        t: Epoch,
        sv: SV,
        rx_m: &Vector3<f64>,
        elev_deg: f64,
        az_deg: f64,
    ) -> f64 {
        let rtm = Self::bias_runtime(t, sv, rx_m, elev_deg, az_deg);
        let modeled = TroposphereModel::Niel.bias_m(&rtm);
        let residual =
            self.inner.scenario.tropo.residual_zwd_m / elev_deg.to_radians().sin().max(0.1);
        modeled + residual
    }

    fn bias_runtime(
        t: Epoch,
        sv: SV,
        rx_m: &Vector3<f64>,
        elev_deg: f64,
        az_deg: f64,
    ) -> BiasRuntime {
        let (lat, lon, h) = ecef_to_geodetic(rx_m);
        BiasRuntime {
            sv,
            epoch: t,
            sv_position_ecef_m: (0.0, 0.0, 0.0),
            sv_elevation_azimuth_deg_deg: (elev_deg, az_deg),
            rx_position_m: *rx_m,
            rcvr_lat_long_alt_deg_deg_km: (lat.to_degrees(), lon.to_degrees(), h / 1.0E3),
            frequency_hz: Carrier::L1.frequency_hz(),
        }
    }

    /// Integer ambiguity (cycles) of a station/satellite/carrier, including cycle slips
    /// that occurred before `t_s` (rover only).
    pub fn ambiguity(&self, station: Station, sv: SV, second: bool, t_s: f64) -> i64 {
        let salt = station.salt() * 7919
            + sv.prn as u64 * 104_729
            + (sv.constellation as u64) * 31
            + second as u64;
        let mut rng = SimRng::new(self.inner.scenario.seed).fork(salt);
        let mut n = rng.uniform_in(-5.0E5, 5.0E5).round() as i64;

        if station == Station::Rover {
            for slip in self.inner.scenario.cycle_slips.iter() {
                if slip.sv == sv && t_s >= slip.at_s {
                    n += if second { slip.cycles_2 } else { slip.cycles_1 };
                }
            }
        }

        n
    }

    fn in_outage(&self, sv: SV, t_s: f64) -> bool {
        self.inner
            .scenario
            .outages
            .iter()
            .any(|o| o.sv == sv && t_s >= o.start_s && t_s < o.end_s)
    }

    /// Computes the signal truth for one station / satellite at sample `k`.
    pub fn signal_truth(&self, station: Station, sv: SV, k: usize) -> Option<SignalTruth> {
        let orbit = self.orbit(sv)?;
        let t = self.inner.scenario.epoch(k);
        let t_tag_s = self.seconds_since_start(t);
        let dt_r = self.station_clock_s(station, k);
        let t_rx_s = t_tag_s - dt_r;

        let rx_m = self.station_position_m(station, t_rx_s);

        // light time iteration (Sagnac included)
        let mut tau = 0.075;
        let mut rho = 0.0;
        let mut sat_rot = Vector3::zeros();

        for _ in 0..8 {
            let sat = orbit.ecef_position_m(t_rx_s - tau);
            let (s, c) = (EARTH_ANGULAR_VEL_RAD * tau).sin_cos();
            sat_rot = Vector3::new(c * sat[0] + s * sat[1], -s * sat[0] + c * sat[1], sat[2]);
            rho = (sat_rot - rx_m).norm();
            tau = rho / SPEED_OF_LIGHT_M_S;
        }

        let (elev_deg, az_deg) = elevation_azimuth_deg(&rx_m, &sat_rot);

        let (r_s, r_r) = (sat_rot.norm(), rx_m.norm());
        let shapiro_m = 2.0 * GM_M3_S2 / SPEED_OF_LIGHT_M_S.powi(2)
            * ((r_s + r_r + rho) / (r_s + r_r - rho)).ln();

        Some(SignalTruth {
            sv,
            geometric_range_m: rho,
            sat_position_m: sat_rot,
            shapiro_m,
            time_of_flight_s: tau,
            elevation_deg: elev_deg,
            sat_clock_s: self.sat_clock_s(sv, t),
            tropo_m: self.tropo_delay_m(t, sv, &rx_m, elev_deg, az_deg),
            iono_l1_m: self.iono_delay_m(elev_deg, Carrier::L1.frequency_hz()),
            ambiguity_1: self.ambiguity(station, sv, false, t_tag_s),
            ambiguity_2: self.ambiguity(station, sv, true, t_tag_s),
        })
    }

    /// Generates the observations of one station at sample `k`.
    pub fn observe(&self, station: Station, k: usize) -> Vec<Candidate> {
        let scenario = &self.inner.scenario;
        let t = scenario.epoch(k);
        let t_s = self.seconds_since_start(t);
        let dt_r = self.station_clock_s(station, k);

        let mut out = Vec::with_capacity(16);

        for orbit in self.inner.orbits.iter() {
            let sv = orbit.sv;

            if self.in_outage(sv, t_s) && station == Station::Rover {
                continue;
            }

            let truth = match self.signal_truth(station, sv, k) {
                Some(truth) => truth,
                None => continue,
            };

            if truth.elevation_deg < scenario.elevation_mask_deg {
                continue;
            }

            let (c1, c2) = self.inner.constellation_of[&sv];

            let common = truth.geometric_range_m
                + truth.shapiro_m
                + SPEED_OF_LIGHT_M_S * (dt_r - truth.sat_clock_s)
                + truth.tropo_m;

            let mut noise_rng = SimRng::new(scenario.seed).fork(
                station.salt() * 1_000_003
                    + k as u64 * 7_919
                    + sv.prn as u64 * 17
                    + sv.constellation as u64,
            );

            let mut observations = Vec::with_capacity(2);

            for (index, carrier) in [c1, c2].iter().enumerate() {
                let f_hz = carrier.frequency_hz();
                let lambda = SPEED_OF_LIGHT_M_S / f_hz;
                let iono = self.iono_delay_m(truth.elevation_deg, f_hz);

                let n = if index == 0 {
                    truth.ambiguity_1
                } else {
                    truth.ambiguity_2
                };

                let code = common
                    + iono
                    + noise_rng.normal(scenario.noise.code_sigma(truth.elevation_deg));
                let phase = common - iono
                    + lambda * n as f64
                    + noise_rng.normal(scenario.noise.phase_sigma(truth.elevation_deg));

                // primary signal always has the best SNR
                let snr = 30.0
                    + 15.0 * truth.elevation_deg.to_radians().sin()
                    + if index == 0 { 5.0 } else { 0.0 };

                observations.push(
                    Observation::pseudo_range(*carrier, code, Some(snr))
                        .with_ambiguous_phase_range_m(phase),
                );
            }

            out.push(Candidate::new(sv, t, observations));
        }

        out
    }

    /// Simulated rover epoch
    pub fn rover_epoch(&self, k: usize) -> SimEpoch {
        let t = self.inner.scenario.epoch(k);
        let t_s = self.seconds_since_start(t);
        let dt_r = self.station_clock_s(Station::Rover, k);

        SimEpoch {
            k,
            epoch: t,
            t_s,
            candidates: self.observe(Station::Rover, k),
            truth: Truth {
                position_m: self.station_position_m(Station::Rover, t_s - dt_r),
                clock_offset_s: dt_r,
            },
        }
    }

    /// Iterates all rover epochs
    pub fn rover_epochs(&self) -> impl Iterator<Item = SimEpoch> + '_ {
        (0..self.inner.scenario.num_epochs()).map(move |k| self.rover_epoch(k))
    }

    /// [OrbitSource] implementation
    pub fn orbit_source(&self) -> SimOrbitSource {
        SimOrbitSource { sim: self.clone() }
    }

    /// [SpacebornBias] implementation
    pub fn spaceborn_biases(&self) -> SimSpacebornBias {
        SimSpacebornBias { sim: self.clone() }
    }

    /// [EnvironmentalBias] implementation
    pub fn environment(&self) -> SimEnvironment {
        SimEnvironment { sim: self.clone() }
    }

    /// [RTKBase] implementation (requires a base in the scenario)
    pub fn rtk_base(&self) -> SimBase {
        assert!(
            self.inner.scenario.base.is_some(),
            "scenario has no base station"
        );
        SimBase { sim: self.clone() }
    }
}

/// [OrbitSource] backed by the simulator
pub struct SimOrbitSource {
    sim: Simulator,
}

impl OrbitSource for SimOrbitSource {
    fn state_at(&self, epoch: Epoch, sv: SV, frame: Frame) -> Option<Orbit> {
        *self.sim.inner.orbit_requests.borrow_mut() += 1;
        let orbit = self.sim.orbit(sv)?;
        let t_s = self.sim.seconds_since_start(epoch);
        let r = orbit.ecef_position_m(t_s);
        Some(Orbit::from_position(
            r[0] / 1.0E3,
            r[1] / 1.0E3,
            r[2] / 1.0E3,
            epoch,
            frame,
        ))
    }
}

/// [SpacebornBias] backed by the simulator
pub struct SimSpacebornBias {
    sim: Simulator,
}

impl SpacebornBias for SimSpacebornBias {
    fn clock_bias(&self, rtm: &BiasRuntime) -> SatelliteClockCorrection {
        // the simulated clocks already contain everything: no relativistic correction needed.
        SatelliteClockCorrection::with_relativistic_correction_s(
            self.sim.sat_clock_s(rtm.sv, rtm.epoch),
        )
    }

    fn group_delay(&self, _: &BiasRuntime) -> Duration {
        Duration::ZERO
    }

    fn mw_bias(&self, _: &BiasRuntime) -> f64 {
        0.0
    }
}

/// [EnvironmentalBias] backed by the simulator: Niel troposphere and
/// a (possibly degraded) ionosphere model.
pub struct SimEnvironment {
    sim: Simulator,
}

impl EnvironmentalBias for SimEnvironment {
    fn troposphere_bias_m(&self, rtm: &BiasRuntime) -> f64 {
        TroposphereModel::Niel.bias_m(rtm)
    }

    fn ionosphere_bias_m(&self, rtm: &BiasRuntime) -> f64 {
        let iono = self.sim.inner.scenario.iono;
        let truth = self
            .sim
            .iono_delay_m(rtm.sv_elevation_azimuth_deg_deg.0, rtm.frequency_hz);
        (1.0 - iono.model_error) * truth
    }
}

/// [RTKBase] backed by the simulator
pub struct SimBase {
    sim: Simulator,
}

impl RTKBase for SimBase {
    fn name(&self) -> String {
        "SIMBASE".to_string()
    }

    fn new_epoch(&mut self, _: Epoch) {}

    fn observe(&self, epoch: Epoch) -> Vec<Candidate> {
        let k = self.sim.index_of(epoch);
        self.sim.observe(Station::Base, k)
    }

    fn reference_position_ecef_m(&self, _: Epoch) -> (f64, f64, f64) {
        let b = self.sim.inner.scenario.base.unwrap();
        (b[0], b[1], b[2])
    }
}
