//! Frozen real data: ESBC00DNK (rover) and MOJN00DNK (base), 2020-06-25.
//!
//! Extracted from nav-solutions/data by `scripts/extract_frozen_data.py`:
//! - raw dual frequency observations (GPS L1/L2, Galileo E1/E5a, BeiDou B1I/B3I),
//! - orbital samples: GRG SP3 (GPS, Galileo, 15'), BeiDou broadcast ephemerides (15'),
//!   interpolated (Lagrange) at the instant requested by the solver,
//! - satellite clocks at the observation epochs: GRG CLK (30s) for GPS and Galileo,
//!   broadcast for BeiDou, plus the relativistic correction and group delays.
//!
//! Two datasets: `dense` (1 h, 30 s) and `day` (23.5 h, 5').
//!
//! Each dataset is replayed through the public `Solver` API, the reports reuse the
//! simulator metrics (`tests::sim::runner`). The truth is the antenna reference
//! point (marker + antenna height from the RINEX header).
use std::{
    collections::{BTreeMap, HashMap},
    str::FromStr,
    sync::OnceLock,
};

use nalgebra::Vector3;

use crate::{
    prelude::{
        BiasRuntime, Candidate, Carrier, Constellation, Duration, Epoch, Frame, Observation, Orbit,
        OrbitSource, RTKBase, SV, SatelliteClockCorrection, SpacebornBias,
    },
    tests::sim::geodesy::{ecef_to_enu_matrix, ecef_to_geodetic},
};

mod replay;
mod suite;

pub use replay::{Replay, replay};

/// Frozen dataset
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Dataset {
    /// 2020-06-25 00:00 - 01:00 GPST, 30 s
    Dense,
    /// 2020-06-25 00:00 - 23:30 GPST, 5'
    Day,
}

impl Dataset {
    fn files(&self) -> (&'static str, &'static str, &'static str, &'static str) {
        match self {
            Self::Dense => (
                include_str!("data/dense/esbc00dnk.obs"),
                include_str!("data/dense/mojn00dnk.obs"),
                include_str!("data/dense/orbits.dat"),
                include_str!("data/dense/clocks.dat"),
            ),
            Self::Day => (
                include_str!("data/day/esbc00dnk.obs"),
                include_str!("data/day/mojn00dnk.obs"),
                include_str!("data/day/orbits.dat"),
                include_str!("data/day/clocks.dat"),
            ),
        }
    }

    /// Loads (once) and returns this dataset
    pub fn load(&self) -> &'static FrozenData {
        static DENSE: OnceLock<FrozenData> = OnceLock::new();
        static DAY: OnceLock<FrozenData> = OnceLock::new();

        let cell = match self {
            Self::Dense => &DENSE,
            Self::Day => &DAY,
        };

        cell.get_or_init(|| {
            let (rover, base, orbits, clocks) = self.files();
            FrozenData {
                rover: StationData::parse(rover),
                base: StationData::parse(base),
                orbits: FrozenOrbits::parse(orbits),
                clocks: FrozenClocks::parse(clocks),
            }
        })
    }
}

/// One loaded dataset
pub struct FrozenData {
    pub rover: StationData,
    pub base: StationData,
    pub orbits: FrozenOrbits,
    pub clocks: FrozenClocks,
}

fn parse_t0(text: &str) -> Epoch {
    let line = text
        .lines()
        .find(|l| l.starts_with("# t0 "))
        .expect("frozen data: missing t0");
    Epoch::from_str(line.trim_start_matches("# t0 ").trim()).expect("frozen data: invalid t0")
}

fn parse_sv(s: &str) -> SV {
    let constellation = match &s[..1] {
        "G" => Constellation::GPS,
        "E" => Constellation::Galileo,
        "C" => Constellation::BeiDou,
        c => panic!("frozen data: unsupported constellation {c}"),
    };
    SV::new(
        constellation,
        s[1..].parse().expect("frozen data: invalid PRN"),
    )
}

fn parse_f64(s: &str) -> Option<f64> {
    if s == "-" { None } else { s.parse().ok() }
}

/// Observations of one station
pub struct StationData {
    pub name: String,
    pub marker_ecef_m: Vector3<f64>,
    /// Antenna height, east, north eccentricities (m)
    pub antenna_hen_m: (f64, f64, f64),
    pub t0: Epoch,
    /// Candidates per epoch
    pub epochs: BTreeMap<Epoch, Vec<Candidate>>,
}

impl StationData {
    fn parse(text: &str) -> Self {
        let t0 = parse_t0(text);

        let mut name = String::new();
        let mut marker = Vector3::zeros();
        let mut hen = (0.0, 0.0, 0.0);
        let mut epochs: BTreeMap<Epoch, Vec<Candidate>> = BTreeMap::new();

        for line in text.lines() {
            let items = line.split_whitespace().collect::<Vec<_>>();

            if line.starts_with("# station ") {
                name = items[2].to_string();
            } else if line.starts_with("# marker_ecef_m ") {
                marker = Vector3::new(
                    items[2].parse().unwrap(),
                    items[3].parse().unwrap(),
                    items[4].parse().unwrap(),
                );
            } else if line.starts_with("# antenna_hen_m ") {
                hen = (
                    items[2].parse().unwrap(),
                    items[3].parse().unwrap(),
                    items[4].parse().unwrap(),
                );
            } else if !line.starts_with('#') && !items.is_empty() {
                let epoch = t0 + Duration::from_seconds(items[0].parse().unwrap());
                let sv = parse_sv(items[1]);

                let observations = items[2..]
                    .chunks(4)
                    .map(|sig| {
                        let carrier =
                            Carrier::from_str(sig[0]).expect("frozen data: invalid carrier");
                        let mut obs = Observation::pseudo_range(
                            carrier,
                            sig[1].parse().unwrap(),
                            parse_f64(sig[3]),
                        );
                        if let Some(phase_m) = parse_f64(sig[2]) {
                            obs = obs.with_ambiguous_phase_range_m(phase_m);
                        }
                        obs
                    })
                    .collect();

                epochs
                    .entry(epoch)
                    .or_default()
                    .push(Candidate::new(sv, epoch, observations));
            }
        }

        Self {
            name,
            marker_ecef_m: marker,
            antenna_hen_m: hen,
            t0,
            epochs,
        }
    }

    /// Antenna reference point (ECEF m)
    pub fn arp_ecef_m(&self) -> Vector3<f64> {
        let (lat, lon, _) = ecef_to_geodetic(&self.marker_ecef_m);
        let enu_to_ecef = ecef_to_enu_matrix(lat, lon).transpose();
        let (h, e, n) = self.antenna_hen_m;
        self.marker_ecef_m + enu_to_ecef * Vector3::new(e, n, h)
    }

    /// Candidates of this epoch, restricted to the given constellations
    pub fn candidates(&self, epoch: Epoch, constellations: &[Constellation]) -> Vec<Candidate> {
        self.epochs
            .get(&epoch)
            .map(|cds| {
                cds.iter()
                    .filter(|cd| constellations.contains(&cd.sv.constellation))
                    .cloned()
                    .collect()
            })
            .unwrap_or_default()
    }
}

/// Orbital samples, interpolated on request
pub struct FrozenOrbits {
    samples: HashMap<SV, Vec<(Epoch, [f64; 3])>>,
}

impl FrozenOrbits {
    /// Lagrange interpolation order (number of samples - 1)
    const ORDER: usize = 10;

    fn parse(text: &str) -> Self {
        let t0 = parse_t0(text);
        let mut samples: HashMap<SV, Vec<(Epoch, [f64; 3])>> = HashMap::new();

        for line in text.lines().filter(|l| !l.starts_with('#')) {
            let items = line.split_whitespace().collect::<Vec<_>>();
            if items.len() != 5 {
                continue;
            }
            let epoch = t0 + Duration::from_seconds(items[0].parse().unwrap());
            samples.entry(parse_sv(items[1])).or_default().push((
                epoch,
                [
                    items[2].parse().unwrap(),
                    items[3].parse().unwrap(),
                    items[4].parse().unwrap(),
                ],
            ));
        }

        for list in samples.values_mut() {
            list.sort_by(|a, b| a.0.cmp(&b.0));
        }

        Self { samples }
    }

    /// Interpolated ECEF position (km)
    pub fn position_km(&self, sv: SV, epoch: Epoch) -> Option<[f64; 3]> {
        let list = self.samples.get(&sv)?;

        if list.len() < Self::ORDER + 1 {
            return None;
        }

        // no extrapolation beyond one sampling period
        let (first, last) = (list[0].0, list[list.len() - 1].0);
        let max_gap = Duration::from_seconds(900.0);

        if epoch < first - max_gap || epoch > last + max_gap {
            return None;
        }

        // closest samples
        let pos = list.partition_point(|(t, _)| *t < epoch);
        let start = pos
            .saturating_sub((Self::ORDER + 1) / 2)
            .min(list.len() - (Self::ORDER + 1));

        let window = &list[start..start + Self::ORDER + 1];

        let samples = window
            .iter()
            .map(|(t, p)| ((*t - epoch).to_seconds(), *p))
            .collect::<Vec<_>>();

        let mut out = [0.0; 3];

        for (i, (t_i, p_i)) in samples.iter().enumerate() {
            let mut l_i = 1.0;
            for (j, (t_j, _)) in samples.iter().enumerate() {
                if i != j {
                    l_i *= -t_j / (t_i - t_j);
                }
            }
            for k in 0..3 {
                out[k] += l_i * p_i[k];
            }
        }

        Some(out)
    }
}

impl OrbitSource for &FrozenOrbits {
    fn state_at(&self, epoch: Epoch, sv: SV, frame: Frame) -> Option<Orbit> {
        let p = self.position_km(sv, epoch)?;
        Some(Orbit::from_position(p[0], p[1], p[2], epoch, frame))
    }
}

/// Satellite clocks at the observation epochs
pub struct FrozenClocks {
    /// (clock, relativistic correction, group delay) in seconds
    clocks: HashMap<(SV, Epoch), (f64, f64, f64)>,
}

impl FrozenClocks {
    fn parse(text: &str) -> Self {
        let t0 = parse_t0(text);
        let mut clocks = HashMap::new();

        for line in text.lines().filter(|l| !l.starts_with('#')) {
            let items = line.split_whitespace().collect::<Vec<_>>();
            if items.len() != 5 {
                continue;
            }
            let epoch = t0 + Duration::from_seconds(items[0].parse().unwrap());
            clocks.insert(
                (parse_sv(items[1]), epoch),
                (
                    items[2].parse().unwrap(),
                    items[3].parse().unwrap(),
                    items[4].parse().unwrap(),
                ),
            );
        }

        Self { clocks }
    }

    pub fn get(&self, sv: SV, epoch: Epoch) -> Option<(f64, f64, f64)> {
        self.clocks.get(&(sv, epoch)).copied()
    }
}

impl SpacebornBias for &FrozenClocks {
    fn clock_bias(&self, rtm: &BiasRuntime) -> SatelliteClockCorrection {
        // relativistic correction included here: the solver does not apply it.
        let (clock_s, rel_s, _) = self.get(rtm.sv, rtm.epoch).unwrap_or_default();
        SatelliteClockCorrection::with_relativistic_correction(Duration::from_seconds(
            clock_s + rel_s,
        ))
    }

    fn group_delay(&self, rtm: &BiasRuntime) -> Duration {
        let (_, _, gd_s) = self.get(rtm.sv, rtm.epoch).unwrap_or_default();
        Duration::from_seconds(gd_s)
    }

    fn mw_bias(&self, _: &BiasRuntime) -> f64 {
        0.0
    }
}

/// Reference station
pub struct FrozenBase {
    pub station: &'static StationData,
    pub constellations: Vec<Constellation>,
}

impl RTKBase for FrozenBase {
    fn name(&self) -> String {
        self.station.name.clone()
    }

    fn new_epoch(&mut self, _: Epoch) {}

    fn observe(&self, epoch: Epoch) -> Vec<Candidate> {
        self.station.candidates(epoch, &self.constellations)
    }

    fn reference_position_ecef_m(&self, _: Epoch) -> (f64, f64, f64) {
        let p = self.station.arp_ecef_m();
        (p[0], p[1], p[2])
    }
}
