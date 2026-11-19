use log::{debug, error};

#[cfg(doc)]
use crate::prelude::TimeScale;

mod dop;
mod kalman;
mod postfit;
mod ppp_ar;

pub(crate) mod apriori;
pub(crate) mod solutions;
pub(crate) mod state;
pub(crate) mod sv;
pub(crate) mod vector;

use nalgebra::{DMatrix, DVector, DimName, U3, U4};

use crate::{
    candidate::differences::Differences,
    constants::SPEED_OF_LIGHT_M_S,
    navigation::{
        apriori::Apriori, dop::DilutionOfPrecision, kalman::Kalman, postfit::PostfitKf,
        ppp_ar::ARSolver, state::State, sv::SVContribution,
    },
    prelude::{Candidate, Config, Duration, Epoch, Error, Frame, Method, UserParameters, Vector3},
};

/// Everything the navigation filter needs to form exact double differences
/// in RTK mode.
pub(crate) struct RtkContext<'a> {
    /// Reference (base) station position, ECEF (m)
    pub base_position_ecef_m: Vector3<f64>,

    /// Rover pivot satellite
    pub rover_pivot: &'a Candidate,

    /// Base pivot satellite
    pub base_pivot: &'a Candidate,

    /// Base station candidates (pivot excluded)
    pub base_candidates: &'a [Candidate],

    /// Observed double [Differences]
    pub double_differences: &'a Differences,
}

/// One linearized measurement
struct Row {
    /// candidate index
    index: usize,

    /// observed - computed, at the nominal state, excluding the receiver clock (m)
    residual_m: f64,

    /// partial derivatives with respect to the position (ECEF)
    h: (f64, f64, f64),

    /// measurement variance (m²)
    variance_m2: f64,

    /// [SVContribution]
    contribution: SVContribution,
}

/// Zenith measurement standard deviations (m), for each navigation [Method].
/// The variance model is σ² = σ0² (1 + 1/sin²(elevation)), as in most GNSS
/// processing engines.
fn zenith_sigma_m(method: Method) -> f64 {
    match method {
        // single frequency code, including residual ionosphere & broadcast errors
        Method::SPP => 0.5,
        // iono free code combination (~x3 noise amplification)
        Method::CPP => 1.0,
        // iono free phase combination
        Method::PPP => 0.01,
    }
}

/// Measurement variance (m²) for this elevation
fn measurement_variance_m2(method: Method, elevation_deg: Option<f64>, uses_rtk: bool) -> f64 {
    let sigma_0 = zenith_sigma_m(method);
    let elev_rad = elevation_deg.unwrap_or(90.0).max(5.0).to_radians();
    let mut var = sigma_0.powi(2) * (1.0 + 1.0 / elev_rad.sin().powi(2));

    if uses_rtk {
        // double differences combine four measurements
        var *= 4.0;
    }

    var
}

fn orbit_position_m(candidate: &Candidate) -> Option<Vector3<f64>> {
    let orbit = candidate.orbit?;
    let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;
    Some(Vector3::new(pos_vel_m[0], pos_vel_m[1], pos_vel_m[2]))
}

/// [Navigation] Solver.
///
/// The navigation filter is an extended Kalman filter, expressed in error-state
/// form: the nominal [State] is held outside of the filter, the filter only
/// estimates the correction δx (and its covariance) at each epoch. Once the
/// correction has been applied, the error-state is reset to zero.
///
/// - Absolute (PPP technique): δx = (δx, δy, δz, c.δt)
/// - Differential (RTK technique): δx = (δx, δy, δz)
pub(crate) struct Navigation {
    /// [Config] preset
    cfg: Config,

    /// [Frame]
    frame: Frame,

    /// F
    f_k: DMatrix<f64>,

    /// Q
    q_k: DMatrix<f64>,

    /// Prefit
    ar_prefit: Option<ARSolver>,

    /// [Kalman]
    kalman: Kalman,

    /// Postfit (Kalman)
    postfit: Option<PostfitKf>,

    /// [SVContribution]s
    pub sv: Vec<SVContribution>,

    /// Published [State] (possibly denoised by the postfit filter)
    pub state: State,

    /// Navigation filter nominal [State]. It differs from the published
    /// [State] when postfit denoising is active: the denoiser must not feed
    /// back into the navigation filter linearization point.
    pub nominal: State,

    /// [DilutionOfPrecision]
    pub dop: DilutionOfPrecision,

    /// Null on first iter
    prev_epoch: Option<Epoch>,
}

impl Navigation {
    /// Creates new [Navigation] solver.
    pub fn new(cfg: &Config, frame: Frame) -> Self {
        Self {
            frame,
            postfit: None,
            prev_epoch: None,
            cfg: cfg.clone(),
            ar_prefit: None,
            state: Default::default(),
            nominal: Default::default(),
            sv: Vec::with_capacity(8),
            kalman: Kalman::new(U4::DIM),
            dop: DilutionOfPrecision::default(),
            f_k: DMatrix::<f64>::identity(U4::DIM, U4::DIM),
            q_k: DMatrix::<f64>::zeros(U4::DIM, U4::DIM),
        }
    }

    /// Reset [Navigation] filter.
    ///
    /// After reset, the [Navigation] solver is in the same state as deployment,
    /// ready to consume a first measurement.
    pub fn reset(&mut self) {
        self.sv.clear();
        self.kalman.reset();

        if let Some(prefit) = &mut self.ar_prefit {
            prefit.reset();
        }

        if let Some(postfit) = &mut self.postfit {
            postfit.reset();
        }

        self.prev_epoch = None;
        self.dop = DilutionOfPrecision::default();
    }

    /// Returns clock index
    pub(crate) fn clock_index() -> usize {
        U4::DIM - 1
    }

    /// Mutable iteration of the [Navigation] filter.
    ///
    /// ## Input
    /// - epoch: sampling [Epoch]
    /// - params: [UserParameters]
    /// - initial_state: nominal [State] to linearize about
    /// - candidates: proposed [Candidate]s
    /// - size: number of proposed [Candidate]s
    /// - rtk: [RtkContext] when using the RTK technique
    pub fn solve(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        initial_state: &State,
        candidates: &[Candidate],
        size: usize,
        rtk: Option<RtkContext>,
    ) -> Result<(), Error> {
        self.sv.clear();

        let uses_rtk = rtk.is_some();
        let ndf = if uses_rtk { U3::DIM } else { U4::DIM };

        let dt = match self.prev_epoch {
            Some(past_epoch) => epoch - past_epoch,
            None => Duration::ZERO,
        };

        if !uses_rtk && self.cfg.method == Method::PPP && self.ar_prefit.is_none() {
            self.ar_prefit = Some(ARSolver::new(self.cfg.clone(), self.frame));
        }

        if !uses_rtk {
            if let Some(prefit) = &mut self.ar_prefit {
                match prefit.run(epoch, params, initial_state, candidates, size) {
                    Ok(_) => {},
                    Err(e) => {
                        error!("{epoch} - ppp prefit failed with {e}");
                        return Err(Error::FloatAmbiguitiesSolving);
                    },
                }
            }
        }

        if !self.kalman.initialized {
            self.kf_initialization(epoch, initial_state, candidates, size, &rtk)?;
        } else {
            // context switch (RTK <=> PPP) is allowed at any epoch
            self.kalman_dimension_mut(ndf);

            self.f_k = DMatrix::identity(ndf, ndf);
            self.q_k = DMatrix::zeros(ndf, ndf);
            params.q_matrix(&mut self.q_k, dt, ndf);

            self.kf_run(epoch, initial_state, candidates, size, &rtk)?;
        }

        self.nominal = self.state.clone();

        if self.cfg.solver.postfit_denoising > 0.0 {
            if let Some(postfit) = &mut self.postfit {
                if let Some(prev_epoch) = self.prev_epoch {
                    let dt = epoch - prev_epoch;
                    let dx = postfit.run(&self.state, dt)?;

                    self.state
                        .postfit_update_mut(self.frame, &dx.x)
                        .map_err(|e| {
                            error!(
                                "{epoch} - postfit state update failed with physical error: {e}"
                            );
                            Error::StateUpdate
                        })?;
                }
            } else {
                self.postfit = Some(PostfitKf::new(
                    &self.state,
                    1.0 / self.cfg.solver.postfit_denoising,
                    1.0 / self.cfg.solver.postfit_denoising,
                    1.0,
                    1.0,
                ));
            }
        }

        self.prev_epoch = Some(epoch);

        Ok(())
    }

    /// Adapts the filter dimension, in case of context switch.
    fn kalman_dimension_mut(&mut self, ndf: usize) {
        let current = self.kalman.covariance().nrows();

        if current == ndf {
            return;
        }

        debug!("navigation filter dimension: {current} -> {ndf}");

        let mut p = DMatrix::<f64>::zeros(ndf, ndf);

        for i in 0..ndf.min(current) {
            for j in 0..ndf.min(current) {
                p[(i, j)] = self.kalman.covariance()[(i, j)];
            }
        }

        // newly introduced states are unknown
        for i in current..ndf {
            p[(i, i)] = (1.0E-3 * SPEED_OF_LIGHT_M_S).powi(2);
        }

        self.kalman.initialize_error_state(&p);
    }

    /// Linearizes one measurement about `position_m`.
    fn measurement_row(
        &self,
        t: Epoch,
        index: usize,
        candidate: &Candidate,
        position_m: Vector3<f64>,
        rtk: &Option<RtkContext>,
    ) -> Result<Row, Error> {
        let mut contribution = SVContribution::default();
        contribution.sv = candidate.sv;

        let variance_m2 =
            measurement_variance_m2(self.cfg.method, candidate.elevation_deg, rtk.is_some());

        match rtk {
            Some(rtk) => {
                let observed = candidate
                    .rtk_vector_contribution(
                        t,
                        &self.cfg,
                        rtk.double_differences,
                        &mut contribution,
                    )?
                    .row_1;

                if let Some((elev, azim)) = candidate.attitude() {
                    contribution.elevation_deg = elev;
                    contribution.azimuth_deg = azim;
                }

                let base = rtk
                    .base_candidates
                    .iter()
                    .find(|cd| cd.sv == candidate.sv)
                    .ok_or(Error::RtkDDPostfitMissing)?;

                let r_i = orbit_position_m(candidate).ok_or(Error::UnresolvedState)?;
                let r_p = orbit_position_m(rtk.rover_pivot).ok_or(Error::UnresolvedState)?;
                let b_i = orbit_position_m(base).ok_or(Error::UnresolvedState)?;
                let b_p = orbit_position_m(rtk.base_pivot).ok_or(Error::UnresolvedState)?;

                let x_b = rtk.base_position_ecef_m;

                // exact double difference of geometric ranges
                let computed = ((r_i - position_m).norm() - (r_p - position_m).norm())
                    - ((b_i - x_b).norm() - (b_p - x_b).norm());

                let h = candidate.rtk_matrix_contribution(position_m, (r_p[0], r_p[1], r_p[2]));

                Ok(Row {
                    index,
                    residual_m: observed - computed,
                    h,
                    variance_m2,
                    contribution,
                })
            },
            None => {
                let vec = candidate.ppp_vector_contribution(
                    &self.cfg,
                    false,
                    position_m,
                    &mut contribution,
                )?;

                let mut residual_m = vec.row_1;

                if self.cfg.method == Method::PPP {
                    // phase navigation requires the float ambiguities
                    let prefit = self
                        .ar_prefit
                        .as_ref()
                        .ok_or(Error::FloatAmbiguitiesSolving)?;

                    let ambiguity_m = prefit
                        .float_ambiguity_m(candidate.sv)
                        .ok_or(Error::UnresolvedAmbiguity)?;

                    residual_m -= ambiguity_m;
                }

                let h = candidate.ppp_matrix_contribution(&self.cfg, position_m);

                Ok(Row {
                    index,
                    residual_m,
                    h,
                    variance_m2,
                    contribution,
                })
            },
        }
    }

    /// Linearizes all measurements about `state`.
    fn measurement_rows(
        &self,
        t: Epoch,
        state: &State,
        candidates: &[Candidate],
        size: usize,
        rtk: &Option<RtkContext>,
    ) -> Vec<Row> {
        let position_m = state.to_position_ecef_m();
        let mut rows = Vec::with_capacity(size);

        for (i, candidate) in candidates.iter().take(size).enumerate() {
            match self.measurement_row(t, i, candidate, position_m, rtk) {
                Ok(row) => rows.push(row),
                Err(e) => {
                    error!("{}({}) - measurement error: {}", t, candidate.sv, e);
                },
            }
        }

        rows
    }

    /// Forms (G, W, z) at the nominal `state`.
    /// z = observed - computed (receiver clock included in absolute mode).
    fn linear_system(
        state: &State,
        rows: &[Row],
        uses_rtk: bool,
    ) -> (DMatrix<f64>, DMatrix<f64>, DVector<f64>) {
        let ndf = if uses_rtk { U3::DIM } else { U4::DIM };
        let n = rows.len();

        let mut g = DMatrix::<f64>::zeros(n, ndf);
        let mut w = DMatrix::<f64>::zeros(n, n);
        let mut z = DVector::<f64>::zeros(n);

        let clock_m = if uses_rtk {
            0.0
        } else {
            state.clock_profile_s().0 * SPEED_OF_LIGHT_M_S
        };

        for (i, row) in rows.iter().enumerate() {
            g[(i, 0)] = row.h.0;
            g[(i, 1)] = row.h.1;
            g[(i, 2)] = row.h.2;

            if !uses_rtk {
                g[(i, Self::clock_index())] = 1.0;
            }

            w[(i, i)] = 1.0 / row.variance_m2;
            z[i] = row.residual_m - clock_m;
        }

        (g, w, z)
    }

    /// Applies a correction δx to the nominal `state`.
    fn apply_correction(
        &self,
        t: Epoch,
        state: &mut State,
        dx: &DVector<f64>,
        uses_rtk: bool,
    ) -> Result<(), Error> {
        state
            .spatial_correction_mut(self.frame, (dx[0], dx[1], dx[2]))
            .map_err(|e| {
                error!("{t} - state update failed with physical error: {e}");
                Error::StateUpdate
            })?;

        if !uses_rtk {
            state.clock_correction_mut(dx[Self::clock_index()] / SPEED_OF_LIGHT_M_S);
        }

        Ok(())
    }

    /// Filter first iteration: iterated (Gauss-Newton) weighted least squares.
    pub fn kf_initialization(
        &mut self,
        t: Epoch,
        state: &State,
        candidates: &[Candidate],
        size: usize,
        rtk: &Option<RtkContext>,
    ) -> Result<(), Error> {
        const MAX_ITER: usize = 10;
        const CONVERGENCE_M: f64 = 1.0E-4;

        let uses_rtk = rtk.is_some();
        let ndf = if uses_rtk { U3::DIM } else { U4::DIM };

        let mut pending = state.clone();
        pending.resize_mut(ndf);

        let mut rows = Vec::new();
        let mut p_0 = DMatrix::<f64>::zeros(ndf, ndf);
        let mut gt_g_inv = DMatrix::<f64>::zeros(ndf, ndf);

        for ith in 0..MAX_ITER {
            rows = self.measurement_rows(t, &pending, candidates, size, rtk);

            if rows.len() < ndf {
                return Err(Error::MatrixMinimalDimension);
            }

            let (g, w, z) = Self::linear_system(&pending, &rows, uses_rtk);

            let gt = g.transpose();
            let gt_w_g = &gt * &w * &g;

            p_0 = gt_w_g.try_inverse().ok_or(Error::MatrixInversion)?;
            gt_g_inv = (&gt * &g).try_inverse().ok_or(Error::MatrixInversion)?;

            let dx = &p_0 * &gt * &w * z;

            debug!("(i={ith}) {t} - dx={dx}");

            self.apply_correction(t, &mut pending, &dx, uses_rtk)?;

            let spatial_correction = (dx[0].powi(2) + dx[1].powi(2) + dx[2].powi(2)).sqrt();

            if spatial_correction < CONVERGENCE_M {
                break;
            }
        }

        let dop = DilutionOfPrecision::new(&pending, gt_g_inv);
        self.state_validation(&dop)?;

        self.kalman.initialize_error_state(&p_0);

        self.sv = rows.into_iter().map(|row| row.contribution).collect();
        self.state = pending;
        self.dop = dop;

        debug!("{} - new state {}", t, self.state);
        debug!("{} - gdop={} tdop={}", t, self.dop.gdop, self.dop.tdop);

        Ok(())
    }

    /// [Kalman] filter run (error-state EKF).
    pub fn kf_run(
        &mut self,
        t: Epoch,
        nominal: &State,
        candidates: &[Candidate],
        size: usize,
        rtk: &Option<RtkContext>,
    ) -> Result<(), Error> {
        let uses_rtk = rtk.is_some();
        let ndf = if uses_rtk { U3::DIM } else { U4::DIM };

        let mut pending = nominal.clone();
        pending.resize_mut(ndf);

        let rows = self.measurement_rows(t, &pending, candidates, size, rtk);

        if rows.len() < ndf {
            return Err(Error::MatrixMinimalDimension);
        }

        let (g, w, z) = Self::linear_system(&pending, &rows, uses_rtk);

        debug!("{t} - innovation: {z}");

        // P⁻ = F P Fᵀ + Q
        self.kalman.predict_error_state(&self.f_k, &self.q_k);

        let estimate = self.kalman.update_error_state(&g, &w, &z)?;

        debug!("{t} - state correction: dx={}", estimate.x);

        self.apply_correction(t, &mut pending, &estimate.x, uses_rtk)?;

        let gt_g_inv = (g.transpose() * &g)
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        let dop = DilutionOfPrecision::new(&pending, gt_g_inv);

        self.state_validation(&dop)?;

        self.dop = dop;
        self.state = pending;
        self.sv = rows.into_iter().map(|row| row.contribution).collect();

        debug!("{t} - new state {}", self.state);
        debug!("{} - gdop={} tdop={}", t, self.dop.gdop, self.dop.tdop);

        Ok(())
    }

    /// Validate pending [State]
    fn state_validation(&self, dop: &DilutionOfPrecision) -> Result<(), Error> {
        if dop.gdop > self.cfg.solver.max_gdop {
            return Err(Error::MaxGdopExceeded);
        }
        Ok(())
    }

    /// True if this filter is initialized
    pub fn is_initialized(&self) -> bool {
        self.kalman.initialized
    }
}

#[cfg(test)]
mod test {
    use crate::navigation::{DilutionOfPrecision, Navigation, State, measurement_variance_m2};
    use crate::prelude::Method;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn navigation_dimensions_clock_index() {
        assert_eq!(Navigation::clock_index(), 3, "invalid clock index for U4");
    }

    #[test]
    fn dilution_of_navigation_precision() {
        let state = State::default();
        let matrix = DMatrix::from_diagonal(&DVector::from_row_slice(&[1.0, 2.0, 3.0, 4.0]));
        let dop = DilutionOfPrecision::new(&state, matrix);
        assert_eq!(dop.gdop, (1.0_f64 + 2.0_f64 + 3.0_f64 + 4.0_f64).sqrt());
        assert_eq!(dop.tdop, 2.0);
    }

    #[test]
    fn measurement_variance_model() {
        for method in [Method::SPP, Method::CPP, Method::PPP] {
            let zenith = measurement_variance_m2(method, Some(90.0), false);
            let low = measurement_variance_m2(method, Some(15.0), false);
            let rtk = measurement_variance_m2(method, Some(90.0), true);
            assert!(low > 5.0 * zenith, "low elevation must be down weighted");
            assert!((rtk - 4.0 * zenith).abs() < 1E-12);
            // unknown elevation behaves like zenith
            assert_eq!(measurement_variance_m2(method, None, false), zenith);
        }
        assert!(
            measurement_variance_m2(Method::PPP, Some(30.0), false)
                < measurement_variance_m2(Method::CPP, Some(30.0), false) / 1000.0
        );
    }
}
