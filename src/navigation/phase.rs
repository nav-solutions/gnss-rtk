//! Carrier phase tracking: one continuous phase arc per satellite.
//!
//! A phase arc starts when a satellite provides dual frequency phase (and code)
//! measurements, and lasts as long as the phase is tracked continuously. Its
//! (float, iono free) ambiguity is a state of the navigation filter. An arc ends,
//! and its ambiguity is reset, when:
//! - the satellite is not observed for longer than the gap tolerance
//!   (data gap, masking, the satellite set),
//! - a cycle slip is detected: geometry free phase jump (with respect to its
//!   linear prediction: the ionosphere varies), or Melbourne-Wübbena jump with
//!   respect to its running mean,
//! - the navigation filter rejects the phase measurement (innovation test).
use log::{debug, info};

use crate::prelude::{Candidate, Duration, Epoch, SV};

/// Minimal gap tolerance
const MIN_GAP_TOLERANCE_S: f64 = 120.0;

/// Gap tolerance, in number of sampling periods
const GAP_TOLERANCE_PERIODS: f64 = 4.0;

/// Geometry free threshold at zenith (m), scaled by 1/sin(elevation): the
/// phase noise (and multipath) grows at low elevation. A (1, 1) cycles slip
/// (GF jump of 5.4 cm on GPS L1/L2) is detected above ~55° (30 s sampling); below, the slip is
/// caught by the navigation filter innovation test (iono free jump of 11 cm).
const GF_THRESHOLD_ZENITH_M: f64 = 0.035;

/// Lowest elevation considered for the GF threshold scaling (deg)
const GF_THRESHOLD_MIN_ELEVATION_DEG: f64 = 10.0;

/// Geometry free threshold (m/s): ionosphere variation allowance, when the
/// ionosphere rate is known (GF predicted from the two previous epochs)
const GF_RATE_THRESHOLD_M_S: f64 = 1.0E-4;

/// Geometry free threshold (m/s²): ionosphere acceleration allowance (long
/// sampling periods: 18 cm over 5', 2 mm over 30 s)
const GF_ACCEL_THRESHOLD_M_S2: f64 = 2.0E-6;

/// Geometry free threshold (m/s): ionosphere drift allowance, second epoch of
/// an arc (ionosphere rate unknown yet)
const GF_FIRST_RATE_THRESHOLD_M_S: f64 = 1.0E-3;

/// Melbourne-Wübbena threshold (m) with respect to the running mean
const MW_THRESHOLD_M: f64 = 3.0;

/// Phase measurements extracted from a [Candidate]
#[derive(Debug, Clone, Copy)]
pub(crate) struct PhaseObservables {
    /// Iono free phase (m)
    pub phase_if_m: f64,
    /// Iono free code (m)
    pub code_if_m: f64,
    /// Geometry free phase (m)
    pub gf_m: f64,
    /// Melbourne-Wübbena (m)
    pub mw_m: f64,
    /// Elevation (deg)
    pub elevation_deg: f64,
}

impl PhaseObservables {
    pub fn from_candidate(cd: &Candidate) -> Option<Self> {
        let phase_if = cd.phase_if_combination()?;
        let code_if = cd.code_if_combination()?;
        let gf = cd.phase_gf_combination()?;
        let wl = cd.phase_wl_combination()?;
        let nl = cd.code_nl_combination()?;

        Some(Self {
            phase_if_m: phase_if.value,
            code_if_m: code_if.value,
            gf_m: gf.value,
            mw_m: wl.value - nl.value,
            elevation_deg: cd.elevation_deg.unwrap_or(90.0),
        })
    }
}

/// Continuous phase arc
#[derive(Debug, Clone)]
pub(crate) struct PhaseArc {
    pub sv: SV,
    /// First epoch of this arc
    pub start: Epoch,
    /// Last epoch this satellite was tracked
    pub last: Epoch,
    /// Number of epochs in this arc
    pub epochs: usize,
    /// Nominal float iono free ambiguity (m)
    pub ambiguity_m: f64,
    /// Last geometry free value (m)
    gf_m: f64,
    /// Geometry free rate (m/s), ionosphere variation
    gf_rate_m_s: Option<f64>,
    /// Melbourne-Wübbena running mean (m)
    mw_mean_m: f64,
}

impl PhaseArc {
    fn new(sv: SV, t: Epoch, obs: &PhaseObservables) -> Self {
        Self {
            sv,
            start: t,
            last: t,
            epochs: 1,
            // initial guess: the code gives the ambiguity at the code noise level
            ambiguity_m: obs.phase_if_m - obs.code_if_m,
            gf_m: obs.gf_m,
            gf_rate_m_s: None,
            mw_mean_m: obs.mw_m,
        }
    }
}

/// Outcome of a tracker update: how to rebuild the filter ambiguity states.
#[derive(Debug, Clone, Default)]
pub(crate) struct TrackerUpdate {
    /// For each arc of the new list: its index in the previous list (None: new arc)
    pub previous_index: Vec<Option<usize>>,
}

/// Phase tracker
#[derive(Debug, Clone, Default)]
pub(crate) struct PhaseTracker {
    /// Arcs, in filter state order
    pub arcs: Vec<PhaseArc>,
    /// Latest epoch
    last_epoch: Option<Epoch>,
    /// Latest sampling period
    sampling: Option<Duration>,
}

impl PhaseTracker {
    pub fn reset(&mut self) {
        self.arcs.clear();
        self.last_epoch = None;
        self.sampling = None;
    }

    pub fn len(&self) -> usize {
        self.arcs.len()
    }

    /// Arc index of this satellite
    pub fn index_of(&self, sv: SV) -> Option<usize> {
        self.arcs.iter().position(|arc| arc.sv == sv)
    }

    fn gap_tolerance(&self) -> Duration {
        let periods = self
            .sampling
            .map(|dt| dt.to_seconds() * GAP_TOLERANCE_PERIODS)
            .unwrap_or(0.0);

        Duration::from_seconds(periods.max(MIN_GAP_TOLERANCE_S))
    }

    /// New epoch: updates the arcs with the current candidates.
    pub fn update(&mut self, t: Epoch, candidates: &[Candidate]) -> TrackerUpdate {
        if let Some(last) = self.last_epoch.filter(|last| t > *last) {
            self.sampling = Some(t - last);
        }

        self.last_epoch = Some(t);

        let gap_tolerance = self.gap_tolerance();

        let mut arcs: Vec<PhaseArc> = Vec::with_capacity(self.arcs.len() + 4);
        let mut previous_index: Vec<Option<usize>> = Vec::with_capacity(self.arcs.len() + 4);

        // existing arcs: continue, slip or lost
        for (index, arc) in self.arcs.iter().enumerate() {
            let observed = candidates
                .iter()
                .find(|cd| cd.sv == arc.sv)
                .and_then(PhaseObservables::from_candidate);

            match observed {
                Some(obs) => {
                    let dt_s = (t - arc.last).to_seconds();

                    if dt_s > gap_tolerance.to_seconds() {
                        info!("{t}({}) - phase arc reset (data gap: {dt_s} s)", arc.sv);
                        arcs.push(PhaseArc::new(arc.sv, t, &obs));
                        previous_index.push(None);
                        continue;
                    }

                    // GF prediction: linear ionosphere variation
                    let (gf_predicted_m, gf_threshold) = match arc.gf_rate_m_s {
                        Some(rate) => (
                            arc.gf_m + rate * dt_s,
                            gf_threshold_m(obs.elevation_deg, dt_s, GF_RATE_THRESHOLD_M_S),
                        ),
                        None => (
                            arc.gf_m,
                            gf_threshold_m(obs.elevation_deg, dt_s, GF_FIRST_RATE_THRESHOLD_M_S),
                        ),
                    };

                    let gf_jump = (obs.gf_m - gf_predicted_m).abs();
                    let mw_jump = (obs.mw_m - arc.mw_mean_m).abs();

                    if gf_jump > gf_threshold || mw_jump > MW_THRESHOLD_M {
                        info!(
                            "{t}({}) - cycle slip (gf jump={:.3} m, mw jump={:.3} m)",
                            arc.sv, gf_jump, mw_jump
                        );
                        arcs.push(PhaseArc::new(arc.sv, t, &obs));
                        previous_index.push(None);
                        continue;
                    }

                    let mut arc = arc.clone();
                    arc.epochs += 1;
                    arc.last = t;
                    arc.gf_rate_m_s = Some((obs.gf_m - arc.gf_m) / dt_s);
                    arc.gf_m = obs.gf_m;
                    arc.mw_mean_m += (obs.mw_m - arc.mw_mean_m) / arc.epochs as f64;

                    arcs.push(arc);
                    previous_index.push(Some(index));
                },
                None => {
                    // not tracked at this epoch: keep the ambiguity while within tolerance
                    if (t - arc.last) <= gap_tolerance {
                        arcs.push(arc.clone());
                        previous_index.push(Some(index));
                    } else {
                        debug!(
                            "{t}({}) - phase arc ended ({} epochs since {})",
                            arc.sv, arc.epochs, arc.start
                        );
                    }
                },
            }
        }

        // new arcs (rising satellites, first phase measurements)
        for cd in candidates.iter() {
            if arcs.iter().any(|arc| arc.sv == cd.sv) {
                continue;
            }

            if let Some(obs) = PhaseObservables::from_candidate(cd) {
                debug!("{t}({}) - new phase arc", cd.sv);
                arcs.push(PhaseArc::new(cd.sv, t, &obs));
                previous_index.push(None);
            }
        }

        self.arcs = arcs;

        TrackerUpdate { previous_index }
    }

    /// Restarts the arc of this satellite (e.g. phase measurement rejected
    /// by the navigation filter): new ambiguity from the current measurements.
    pub fn restart(&mut self, t: Epoch, sv: SV, obs: &PhaseObservables) {
        if let Some(index) = self.index_of(sv) {
            self.arcs[index] = PhaseArc::new(sv, t, obs);
        }
    }
}

/// Geometry free cycle slip threshold (m), at this elevation, `dt_s` after
/// the previous measurement, given the ionosphere rate allowance (m/s)
fn gf_threshold_m(elevation_deg: f64, dt_s: f64, rate_allowance_m_s: f64) -> f64 {
    let sin_el = elevation_deg
        .max(GF_THRESHOLD_MIN_ELEVATION_DEG)
        .to_radians()
        .sin();

    GF_THRESHOLD_ZENITH_M / sin_el
        + rate_allowance_m_s * dt_s
        + GF_ACCEL_THRESHOLD_M_S2 * dt_s.powi(2)
}

/// Initial variance (m²) of a new iono free ambiguity: code noise level.
pub(crate) fn initial_ambiguity_variance_m2() -> f64 {
    10.0_f64.powi(2)
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::prelude::{Carrier, Constellation, Observation};

    use std::str::FromStr;

    const L1_HZ: f64 = 1575.42E6;
    const L2_HZ: f64 = 1227.60E6;
    const C: f64 = crate::constants::SPEED_OF_LIGHT_M_S;

    fn sv(prn: u8) -> SV {
        SV::new(Constellation::GPS, prn)
    }

    fn t0() -> Epoch {
        Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap()
    }

    /// GPS L1/L2 candidate: geometric range `rho`, slant ionosphere (L1) `iono`,
    /// phase cycle offsets (n1, n2), at `elevation_deg`.
    fn candidate(
        prn: u8,
        t: Epoch,
        rho: f64,
        iono: f64,
        n1: i64,
        n2: i64,
        elevation_deg: f64,
    ) -> Candidate {
        let (l1, l2) = (C / L1_HZ, C / L2_HZ);
        let iono_2 = iono * (L1_HZ / L2_HZ).powi(2);

        let mut cd = Candidate::new(
            sv(prn),
            t,
            vec![
                Observation::pseudo_range(Carrier::L1, rho + iono, Some(45.0))
                    .with_ambiguous_phase_range_m(rho - iono + n1 as f64 * l1),
                Observation::pseudo_range(Carrier::L2, rho + iono_2, Some(45.0))
                    .with_ambiguous_phase_range_m(rho - iono_2 + n2 as f64 * l2),
            ],
        );

        cd.elevation_deg = Some(elevation_deg);
        cd
    }

    fn at(k: usize) -> Epoch {
        t0() + Duration::from_seconds(30.0 * k as f64)
    }

    /// Runs `epochs` epochs of one satellite, the slip (n1, n2) occurring at `slip_k`.
    fn run_slip(n1: i64, n2: i64, elevation_deg: f64) -> Vec<TrackerUpdate> {
        let mut tracker = PhaseTracker::default();

        (0..10)
            .map(|k| {
                let (s1, s2) = if k >= 5 { (n1, n2) } else { (0, 0) };
                // slowly varying geometry and ionosphere
                let rho = 2.2E7 + 500.0 * k as f64;
                let iono = 3.0 + 0.001 * k as f64;
                tracker.update(
                    at(k),
                    &[candidate(
                        1,
                        at(k),
                        rho,
                        iono,
                        100 + s1,
                        80 + s2,
                        elevation_deg,
                    )],
                )
            })
            .collect()
    }

    fn slipped(updates: &[TrackerUpdate]) -> Vec<usize> {
        updates
            .iter()
            .enumerate()
            .skip(1)
            .filter(|(_, u)| u.previous_index[0].is_none())
            .map(|(k, _)| k)
            .collect()
    }

    #[test]
    fn continuous_arc() {
        let updates = run_slip(0, 0, 45.0);
        assert!(slipped(&updates).is_empty());
        assert!(updates[0].previous_index[0].is_none(), "new arc");
        assert_eq!(updates[9].previous_index, vec![Some(0)]);
    }

    #[test]
    fn single_carrier_slips_are_detected() {
        for (n1, n2) in [(1, 0), (0, 1), (-1, 0), (0, -2), (5, 0)] {
            assert_eq!(
                slipped(&run_slip(n1, n2, 45.0)),
                vec![5],
                "slip ({n1}, {n2})"
            );
        }
    }

    #[test]
    fn wide_lane_neutral_slip_is_detected_by_gf() {
        // (1, 1): GF jump of 5.4 cm, MW unchanged
        assert_eq!(slipped(&run_slip(1, 1, 70.0)), vec![5]);
    }

    #[test]
    fn gf_neutral_slip_is_detected_by_mw() {
        // (9, 7): GF jump of 3 mm only, wide lane jump of 2 cycles (1.7 m)...
        // ... below the MW threshold: (77, 60) is GF neutral too, with a 17 cycles
        // (14.6 m) wide lane jump
        assert_eq!(slipped(&run_slip(77, 60, 45.0)), vec![5]);
    }

    #[test]
    fn fast_ionosphere_is_not_a_slip() {
        // midday ionosphere, 5' sampling: the GF varies by ~20 cm between epochs
        let mut tracker = PhaseTracker::default();
        let t = |k: usize| t0() + Duration::from_seconds(300.0 * k as f64);

        for k in 0..12 {
            let iono = 3.0 + 0.3 * k as f64 + 0.002 * (k * k) as f64;
            let update = tracker.update(t(k), &[candidate(1, t(k), 2.2E7, iono, 10, 8, 30.0)]);
            if k > 0 {
                assert_eq!(update.previous_index, vec![Some(0)], "false slip at {k}");
            }
        }
    }

    #[test]
    fn gf_threshold_grows_at_low_elevation() {
        let r = GF_RATE_THRESHOLD_M_S;
        assert!((gf_threshold_m(90.0, 0.0, r) - GF_THRESHOLD_ZENITH_M).abs() < 1.0E-12);
        assert!(gf_threshold_m(30.0, 0.0, r) > 1.9 * GF_THRESHOLD_ZENITH_M);
        // floor of the elevation scaling
        assert_eq!(gf_threshold_m(2.0, 0.0, r), gf_threshold_m(10.0, 0.0, r));
        // ionosphere drift allowance
        assert!(gf_threshold_m(90.0, 300.0, r) > gf_threshold_m(90.0, 30.0, r));
    }

    #[test]
    fn data_gaps() {
        let mut tracker = PhaseTracker::default();
        let cd = |k: usize| candidate(1, at(k), 2.2E7, 3.0, 100, 80, 45.0);

        for k in 0..4 {
            tracker.update(at(k), &[cd(k)]);
        }

        // short masking (2 epochs): the arc and its ambiguity are kept
        tracker.update(at(4), &[]);
        assert_eq!(tracker.len(), 1);
        tracker.update(at(5), &[]);
        let update = tracker.update(at(6), &[cd(6)]);
        assert_eq!(update.previous_index, vec![Some(0)]);
        assert_eq!(tracker.arcs[0].start, at(0));

        // long gap (> 4 sampling periods, > 120 s): the arc ends...
        for k in 7..13 {
            tracker.update(at(k), &[]);
        }
        assert_eq!(tracker.len(), 0, "arc ended");

        // ... and a new one starts
        let update = tracker.update(at(13), &[cd(13)]);
        assert_eq!(update.previous_index, vec![None]);
        assert_eq!(tracker.arcs[0].start, at(13));
    }

    #[test]
    fn rising_and_setting_satellites_reindexing() {
        let mut tracker = PhaseTracker::default();
        let cd =
            |prn: u8, k: usize| candidate(prn, at(k), 2.2E7 + prn as f64 * 1.0E5, 3.0, 10, 8, 45.0);

        tracker.update(at(0), &[cd(1, 0), cd(2, 0), cd(3, 0)]);
        assert_eq!(tracker.len(), 3);

        // G02 sets (and is not seen anymore), G04 rises
        let mut last = TrackerUpdate::default();
        for k in 1..8 {
            last = tracker.update(at(k), &[cd(1, k), cd(3, k), cd(4, k)]);
        }

        assert_eq!(
            tracker.arcs.iter().map(|arc| arc.sv).collect::<Vec<_>>(),
            vec![sv(1), sv(3), sv(4)]
        );
        // G01 and G03 continue (old indexes 0 and 1 after G02 left), G04 is not new anymore
        assert_eq!(last.previous_index, vec![Some(0), Some(1), Some(2)]);

        assert_eq!(tracker.index_of(sv(3)), Some(1));
        assert_eq!(tracker.index_of(sv(2)), None);
    }

    #[test]
    fn restart_and_initial_ambiguity() {
        let mut tracker = PhaseTracker::default();
        let cd = candidate(1, at(0), 2.2E7, 3.0, 100, 80, 45.0);

        tracker.update(at(0), std::slice::from_ref(&cd));

        // initial ambiguity: iono free phase - iono free code
        let obs = PhaseObservables::from_candidate(&cd).unwrap();
        assert!((tracker.arcs[0].ambiguity_m - (obs.phase_if_m - obs.code_if_m)).abs() < 1.0E-6);

        // the iono free ambiguity of (100, 80) cycles on L1 / L2
        let (l1, l2) = (C / L1_HZ, C / L2_HZ);
        let (f1, f2) = (L1_HZ * L1_HZ, L2_HZ * L2_HZ);
        let expected = (f1 * 100.0 * l1 - f2 * 80.0 * l2) / (f1 - f2);
        assert!((tracker.arcs[0].ambiguity_m - expected).abs() < 1.0E-4);

        tracker.arcs[0].ambiguity_m = 0.0;
        tracker.arcs[0].epochs = 10;
        tracker.restart(at(1), sv(1), &obs);

        assert_eq!(tracker.arcs[0].epochs, 1);
        assert!((tracker.arcs[0].ambiguity_m - expected).abs() < 1.0E-4);
    }
}
