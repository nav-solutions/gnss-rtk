//! Precise Point Positioning (float): one navigation filter, code and phase.
//!
//! State (error-state EKF):
//!
//! ```text
//! [ dx, dy, dz, c.dt, ZWD, ISB_1, ISB_2, ISB_3, N_1 ... N_n ]
//! ```
//!
//! - ZWD: residual zenith wet delay (the a priori troposphere model being applied),
//! - ISB_k: inter-system bias of the k-th secondary constellation, with respect to
//!   the reference constellation (GPS when available) that defines c.dt,
//! - N_i: float iono free ambiguity of the i-th phase arc ([PhaseTracker]).
//!
//! Each epoch processes the iono free code AND phase of every satellite:
//! - at start, the ambiguities are unknown (code noise level): the solution
//!   is a CPP solution,
//! - as the geometry changes, the filter separates position and ambiguities:
//!   the phase takes over progressively, the solution converges to the PPP
//!   accuracy. No mode switch, no second filter.
//!
//! The phase tracker handles rising / setting satellites, data gaps and cycle
//! slips; the filter also rejects (and restarts) phase arcs whose innovation is
//! inconsistent with its covariance.
use log::{debug, error, info};

use nalgebra::{DMatrix, DVector, U4};

use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    navigation::{
        DilutionOfPrecision, Navigation, measurement_variance_m2,
        phase::{PhaseObservables, initial_ambiguity_variance_m2},
        state::State,
        sv::SVContribution,
    },
    prelude::{Candidate, Duration, Epoch, Error, Method, UserParameters, Vector3},
};

use nalgebra::DimName;

/// Receiver clock index
const CLOCK: usize = 3;

/// Zenith wet delay index
const ZWD: usize = 4;

/// First inter-system bias index
const ISB: usize = 5;

/// Number of inter-system bias states (secondary constellations)
const MAX_ISB: usize = 3;

/// First ambiguity index
const AMBIGUITY: usize = ISB + MAX_ISB;

/// Initial inter-system bias standard deviation (m): ~100 ns
const ISB_INITIAL_SIGMA_M: f64 = 30.0;

/// Inter-system bias random walk (m²/s): stable hardware delays
const ISB_RANDOM_WALK_M2_S: f64 = 1.0E-6;

/// Receiver clock (white noise) standard deviation, after re-centering (m)
const RECEIVER_CLOCK_SIGMA_M: f64 = 100.0;

/// Initial ZWD standard deviation (m)
const ZWD_INITIAL_SIGMA_M: f64 = 0.3;

/// ZWD random walk (m²/s): ~2 cm/sqrt(h)
const ZWD_RANDOM_WALK_M2_S: f64 = 1.0E-7;

/// Phase innovation test: number of standard deviations
const PHASE_INNOVATION_SIGMAS: f64 = 5.0;

/// Phase innovation test: minimal rejection threshold (m). Below, a slip
/// is not distinguishable from multipath on real data.
const PHASE_INNOVATION_MIN_M: f64 = 0.08;

/// Wet mapping function (cosecant law, adequate above the default elevation mask)
fn wet_mapping(elevation_deg: f64) -> f64 {
    1.0 / elevation_deg.max(5.0).to_radians().sin()
}

/// Linearized PPP measurements of one satellite
struct PppRow {
    /// partial derivatives with respect to the position
    h: (f64, f64, f64),
    /// wet mapping
    m_wet: f64,
    /// code: observed - computed (receiver clock, ZWD excluded)
    code_residual_m: f64,
    code_variance_m2: f64,
    /// phase: (observed - computed, variance, arc index), ambiguity excluded
    phase: Option<(f64, f64, usize)>,
    /// phase observables (arc restart)
    observables: Option<PhaseObservables>,
    /// inter-system bias: (state index, current value in m), None for the reference
    isb: Option<(usize, f64)>,
    contribution: SVContribution,
}

impl PppRow {
    /// Current inter-system bias (m)
    fn isb_m(&self) -> f64 {
        self.isb.map(|(_, isb_m)| isb_m).unwrap_or(0.0)
    }
}

impl Navigation {
    /// [Config] used to form the code rows
    fn code_config(&self) -> crate::prelude::Config {
        let mut cfg = self.cfg.clone();
        cfg.method = Method::CPP;
        cfg
    }

    /// PPP rows at the nominal position
    fn ppp_rows(
        &self,
        t: Epoch,
        candidates: &[Candidate],
        size: usize,
        position_m: Vector3<f64>,
    ) -> Vec<PppRow> {
        let code_cfg = self.code_config();
        let mut rows = Vec::with_capacity(size);

        for candidate in candidates.iter().take(size) {
            if !self
                .isb_systems
                .iter()
                .any(|(c, _)| *c == candidate.sv.constellation)
            {
                continue; // no clock / inter-system bias state for this constellation
            }

            let mut contribution = SVContribution {
                sv: candidate.sv,
                ..Default::default()
            };

            let code = match candidate.ppp_vector_contribution(
                &code_cfg,
                false,
                position_m,
                &mut contribution,
            ) {
                Ok(vec) => vec.row_1,
                Err(e) => {
                    error!("{t}({}) - code measurement error: {e}", candidate.sv);
                    continue;
                },
            };

            let elevation_deg = candidate.elevation_deg;
            let h = candidate.ppp_matrix_contribution(&code_cfg, position_m);

            let observables = PhaseObservables::from_candidate(candidate);

            let phase = match (self.tracker.index_of(candidate.sv), observables) {
                (Some(arc), Some(_)) => {
                    let mut c = SVContribution::default();
                    match candidate.ppp_vector_contribution(&self.cfg, true, position_m, &mut c) {
                        Ok(vec) => Some((
                            vec.row_2,
                            measurement_variance_m2(Method::PPP, elevation_deg, false),
                            arc,
                        )),
                        Err(_) => None,
                    }
                },
                _ => None,
            };

            let isb = self
                .isb_systems
                .iter()
                .enumerate()
                .skip(1)
                .find(|(_, (c, _))| *c == candidate.sv.constellation)
                .map(|(slot, (_, isb_m))| (ISB + slot - 1, *isb_m));

            rows.push(PppRow {
                isb,
                h,
                m_wet: wet_mapping(elevation_deg.unwrap_or(90.0)),
                code_residual_m: code,
                code_variance_m2: measurement_variance_m2(Method::CPP, elevation_deg, false),
                phase,
                observables,
                contribution,
            });
        }

        rows
    }

    /// Covariance with the ambiguity states re-indexed after a tracker update.
    fn reindexed_covariance(p: &DMatrix<f64>, previous_index: &[Option<usize>]) -> DMatrix<f64> {
        let n = AMBIGUITY + previous_index.len();
        let mut out = DMatrix::<f64>::zeros(n, n);

        // old index of each new state
        let map = (0..n)
            .map(|i| {
                if i < AMBIGUITY {
                    Some(i)
                } else {
                    previous_index[i - AMBIGUITY].map(|k| AMBIGUITY + k)
                }
            })
            .collect::<Vec<_>>();

        for i in 0..n {
            match map[i] {
                Some(oi) => {
                    for j in 0..n {
                        if let Some(oj) = map[j] {
                            out[(i, j)] = p[(oi, oj)];
                        }
                    }
                },
                None => out[(i, i)] = initial_ambiguity_variance_m2(),
            }
        }

        out
    }

    /// Resets one state of the covariance (no correlation, given variance)
    fn reset_state_covariance(p: &mut DMatrix<f64>, index: usize, variance: f64) {
        let n = p.nrows();
        for k in 0..n {
            p[(index, k)] = 0.0;
            p[(k, index)] = 0.0;
        }
        p[(index, index)] = variance;
    }

    /// Assigns the inter-system bias states: the reference constellation
    /// (GPS when available) defines the receiver clock, the other ones get
    /// an ISB state (first come first served, [MAX_ISB] at most).
    fn ppp_assign_isb(&mut self, t: Epoch, candidates: &[Candidate]) {
        use crate::prelude::Constellation;

        if self.isb_systems.is_empty() {
            let reference = if candidates
                .iter()
                .any(|cd| cd.sv.constellation == Constellation::GPS)
            {
                Some(Constellation::GPS)
            } else {
                candidates.first().map(|cd| cd.sv.constellation)
            };

            if let Some(reference) = reference {
                self.isb_systems.push((reference, 0.0));
            }
        }

        for cd in candidates.iter() {
            let constellation = cd.sv.constellation;

            if self.isb_systems.iter().any(|(c, _)| *c == constellation) {
                continue;
            }

            if self.isb_systems.len() > MAX_ISB {
                error!(
                    "{t}({}) - too many constellations: no inter-system bias state",
                    cd.sv
                );
                continue;
            }

            debug!("{t} - ppp: {constellation} inter-system bias state");
            self.isb_systems.push((constellation, 0.0));
        }
    }

    /// Completes a (code based) initialization with the troposphere and
    /// the phase arcs of this epoch.
    pub(crate) fn ppp_initialize(&mut self, t: Epoch, candidates: &[Candidate]) {
        self.tracker.reset();
        self.zwd_m = 0.0;
        self.isb_systems.clear();
        self.ppp_assign_isb(t, candidates);

        let update = self.tracker.update(t, candidates);

        let p_4 = self.kalman.covariance().clone();
        let n = AMBIGUITY + update.previous_index.len();

        let mut p = DMatrix::<f64>::zeros(n, n);

        for i in 0..U4::DIM.min(p_4.nrows()) {
            for j in 0..U4::DIM.min(p_4.ncols()) {
                p[(i, j)] = p_4[(i, j)];
            }
        }

        p[(ZWD, ZWD)] = ZWD_INITIAL_SIGMA_M.powi(2);

        for i in ISB..AMBIGUITY {
            p[(i, i)] = ISB_INITIAL_SIGMA_M.powi(2);
        }

        for i in AMBIGUITY..n {
            p[(i, i)] = initial_ambiguity_variance_m2();
        }

        self.kalman.initialize_error_state(&p);

        info!("{t} - ppp: {} phase arcs initialized", self.tracker.len());
    }

    /// PPP filter run, transactional (the filter, the tracker and the
    /// troposphere are restored if the epoch is rejected).
    pub(crate) fn ppp_run(
        &mut self,
        t: Epoch,
        params: UserParameters,
        dt: Duration,
        nominal: &State,
        candidates: &[Candidate],
        size: usize,
    ) -> Result<(), Error> {
        let snapshot = (
            self.kalman.clone(),
            self.tracker.clone(),
            self.zwd_m,
            self.isb_systems.clone(),
        );

        let ret = self.ppp_run_uncommitted(t, params, dt, nominal, candidates, size);

        if ret.is_err() {
            (self.kalman, self.tracker, self.zwd_m, self.isb_systems) = snapshot;
        }

        ret
    }

    fn ppp_run_uncommitted(
        &mut self,
        t: Epoch,
        params: UserParameters,
        dt: Duration,
        nominal: &State,
        candidates: &[Candidate],
        size: usize,
    ) -> Result<(), Error> {
        // context switch (from RTK, or another method): restart the PPP states
        if self.kalman.covariance().nrows() < AMBIGUITY {
            self.kalman_dimension_mut(U4::DIM);
            self.ppp_initialize(t, &[]);
        }

        let mut pending = nominal.clone();
        pending.resize_mut(U4::DIM);

        self.ppp_assign_isb(t, candidates);

        // phase tracking: arcs continue, slip, end or start
        let update = self.tracker.update(t, candidates);
        let mut p = Self::reindexed_covariance(self.kalman.covariance(), &update.previous_index);
        let n = p.nrows();

        // prediction
        let mut q = DMatrix::<f64>::zeros(n, n);
        params.q_matrix(&mut q, dt, U4::DIM);
        q[(ZWD, ZWD)] = ZWD_RANDOM_WALK_M2_S * dt.to_seconds().abs();

        for i in ISB..AMBIGUITY {
            q[(i, i)] = ISB_RANDOM_WALK_M2_S * dt.to_seconds().abs();
        }

        p += q;

        // The receiver clock is a white noise process: it is re-centered on the code
        // below, then estimated from scratch at every epoch. Real receiver clocks
        // (steered, or simply noisy) do not follow the user clock profile at the
        // decimeter level: a too tight clock prior turns every common mode phase
        // innovation into a vertical (and troposphere) error.
        for i in 0..n {
            p[(CLOCK, i)] = 0.0;
            p[(i, CLOCK)] = 0.0;
        }
        p[(CLOCK, CLOCK)] = RECEIVER_CLOCK_SIGMA_M.powi(2);

        let rows = self.ppp_rows(t, candidates, size, pending.to_position_ecef_m());

        if rows.len() < U4::DIM + 1 {
            return Err(Error::MatrixMinimalDimension);
        }

        // receiver clock re-centering (median code innovation)
        let clock_m = pending.clock_profile_s().0 * SPEED_OF_LIGHT_M_S;

        let mut code_innov = rows
            .iter()
            .map(|r| r.code_residual_m - clock_m - r.m_wet * self.zwd_m - r.isb_m())
            .collect::<Vec<_>>();

        let mut sorted = code_innov.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let shift_m = sorted[sorted.len() / 2];

        pending.clock_correction_mut(shift_m / SPEED_OF_LIGHT_M_S);
        code_innov.iter_mut().for_each(|z| *z -= shift_m);

        let clock_m = clock_m + shift_m;

        // phase innovations
        let mut phase_innov: Vec<Option<(f64, f64, usize)>> = rows
            .iter()
            .map(|row| {
                row.phase.map(|(residual_m, var_m2, arc)| {
                    let z = residual_m
                        - clock_m
                        - row.m_wet * self.zwd_m
                        - row.isb_m()
                        - self.tracker.arcs[arc].ambiguity_m;
                    (z, var_m2, arc)
                })
            })
            .collect();

        // Phase consistency test (cycle slips missed by the tracker, outliers).
        // The receiver clock is only known at the code level at this point: the
        // test applies to the innovations with respect to their median (common mode
        // removed), and to the state uncertainty without the clock.
        let mature = phase_innov
            .iter()
            .flatten()
            .filter(|(_, _, arc)| self.tracker.arcs[*arc].epochs > 2)
            .map(|(z, _, _)| *z)
            .collect::<Vec<_>>();

        if mature.len() >= 3 {
            let mut sorted = mature.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let common_m = sorted[sorted.len() / 2];

            for (k, row) in rows.iter().enumerate() {
                let Some((z, var_m2, arc)) = phase_innov[k] else {
                    continue;
                };

                if self.tracker.arcs[arc].epochs <= 2 {
                    continue;
                }

                let mut hv = DVector::<f64>::zeros(n);
                hv[0] = row.h.0;
                hv[1] = row.h.1;
                hv[2] = row.h.2;
                hv[ZWD] = row.m_wet;
                if let Some((col, _)) = row.isb {
                    hv[col] = 1.0;
                }
                hv[AMBIGUITY + arc] = 1.0;

                let s = (hv.transpose() * &p * &hv)[(0, 0)] + var_m2;
                let threshold = (PHASE_INNOVATION_SIGMAS * s.sqrt()).max(PHASE_INNOVATION_MIN_M);

                if (z - common_m).abs() <= threshold {
                    continue;
                }

                info!(
                    "{t}({}) - phase innovation {:.3} m > {:.3} m: arc restarted",
                    row.contribution.sv,
                    z - common_m,
                    threshold
                );

                let Some(obs) = row.observables else {
                    phase_innov[k] = None;
                    continue;
                };

                self.tracker.restart(t, row.contribution.sv, &obs);
                Self::reset_state_covariance(
                    &mut p,
                    AMBIGUITY + arc,
                    initial_ambiguity_variance_m2(),
                );

                let (residual_m, _, _) = row.phase.unwrap();
                let z = residual_m
                    - clock_m
                    - row.m_wet * self.zwd_m
                    - row.isb_m()
                    - self.tracker.arcs[arc].ambiguity_m;

                phase_innov[k] = Some((z, var_m2, arc));
            }
        }

        // linear system
        let num_phase = phase_innov.iter().filter(|p| p.is_some()).count();
        let m = rows.len() + num_phase;

        let mut g = DMatrix::<f64>::zeros(m, n);
        let mut w = DMatrix::<f64>::zeros(m, m);
        let mut z = DVector::<f64>::zeros(m);

        // position / clock geometry of the code rows (DOP)
        let mut g_code = DMatrix::<f64>::zeros(rows.len(), U4::DIM);

        let mut i = 0;

        for (k, row) in rows.iter().enumerate() {
            for (col, v) in [(0, row.h.0), (1, row.h.1), (2, row.h.2), (CLOCK, 1.0)] {
                g[(i, col)] = v;
                g_code[(k, col)] = v;
            }
            g[(i, ZWD)] = row.m_wet;
            if let Some((col, _)) = row.isb {
                g[(i, col)] = 1.0;
            }
            w[(i, i)] = 1.0 / row.code_variance_m2;
            z[i] = code_innov[k];
            i += 1;

            if let Some((z_phase, var_m2, arc)) = phase_innov[k] {
                g[(i, 0)] = row.h.0;
                g[(i, 1)] = row.h.1;
                g[(i, 2)] = row.h.2;
                g[(i, CLOCK)] = 1.0;
                g[(i, ZWD)] = row.m_wet;
                if let Some((col, _)) = row.isb {
                    g[(i, col)] = 1.0;
                }
                g[(i, AMBIGUITY + arc)] = 1.0;
                w[(i, i)] = 1.0 / var_m2;
                z[i] = z_phase;
                i += 1;
            }
        }

        debug!(
            "{t} - ppp: {} code, {} phase measurements, {} arcs",
            rows.len(),
            num_phase,
            self.tracker.len()
        );

        self.kalman.initialize_error_state(&p);

        let estimate = self.kalman.update_error_state(&g, &w, &z)?;

        if !estimate
            .x
            .iter()
            .chain(estimate.p.iter())
            .all(|v| v.is_finite())
        {
            error!("{t} - non finite ppp state correction");
            return Err(Error::StateUpdate);
        }

        let dx = &estimate.x;

        debug!(
            "{t} - ppp correction: position=({:.3}, {:.3}, {:.3}) m, clock={:.3} m, zwd={:.4} m",
            dx[0], dx[1], dx[2], dx[CLOCK], dx[ZWD]
        );

        pending
            .spatial_correction_mut(self.frame, (dx[0], dx[1], dx[2]))
            .map_err(|e| {
                error!("{t} - state update failed with physical error: {e}");
                Error::StateUpdate
            })?;

        pending.clock_correction_mut(dx[CLOCK] / SPEED_OF_LIGHT_M_S);

        self.zwd_m += dx[ZWD];

        for (slot, (_, isb_m)) in self.isb_systems.iter_mut().enumerate().skip(1) {
            *isb_m += dx[ISB + slot - 1];
        }

        for (j, arc) in self.tracker.arcs.iter_mut().enumerate() {
            arc.ambiguity_m += dx[AMBIGUITY + j];
        }

        let gt_g_inv = (g_code.transpose() * &g_code)
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        let dop = DilutionOfPrecision::new(&pending, gt_g_inv);

        if dop.gdop > self.cfg.solver.max_gdop {
            return Err(Error::MaxGdopExceeded);
        }

        self.dop = dop;
        self.state = pending;
        self.sv = rows.into_iter().map(|row| row.contribution).collect();

        debug!(
            "{t} - ppp sigma: pos=({:.3}, {:.3}, {:.3}) m, clock={:.3} m, zwd={:.4} m, zwd_nominal={:.4} m",
            estimate.p[(0, 0)].sqrt(),
            estimate.p[(1, 1)].sqrt(),
            estimate.p[(2, 2)].sqrt(),
            estimate.p[(CLOCK, CLOCK)].sqrt(),
            estimate.p[(ZWD, ZWD)].sqrt(),
            self.zwd_m
        );

        Ok(())
    }
}
