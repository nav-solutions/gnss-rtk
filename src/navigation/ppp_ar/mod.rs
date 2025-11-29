#[cfg(doc)]
use crate::prelude::TimeScale;

use log::{debug, error};

use nalgebra::{DMatrix, DVector, DimName, U4, U6, U8};

use crate::{
    candidate::differences::Differences,
    navigation::{
        Navigation,
        dop::DilutionOfPrecision,
        kalman::{Kalman, KfEstimate},
        state::State,
        sv::SVContribution,
    },
    prelude::{Candidate, Config, Duration, Epoch, Error, Frame, Method, SV},
    rtk::RTKBase,
    user::UserParameters,
};

use std::collections::HashMap;

mod lambda;
use lambda::LambdaAR;

/// RTK+PPP Solver
pub struct Solver {
    /// [Config] preset
    cfg: Config,

    /// [Frame]
    frame: Frame,

    /// Y allocation
    y_k_vec: Vec<f64>,

    /// W allocation
    w_k_vec: Vec<f64>,

    /// W
    w_k: DMatrix<f64>,

    /// F
    f_k: DMatrix<f64>,

    /// Q
    q_k: DMatrix<f64>,

    /// G
    g_k: DMatrix<f64>,

    /// X
    x_k: DVector<f64>,

    /// P
    p_k: DMatrix<f64>,

    /// contribution allocation
    indexes: Vec<usize>,

    /// SV contributions
    sv: Vec<SVContribution>,

    /// [Kalman]
    kalman: Kalman,

    /// Lambda X
    lambda_x: DMatrix<f64>,

    /// Lambda Q
    lambda_q: DMatrix<f64>,

    /// Fixed ambiguities
    pub fixed_amb: HashMap<SV, i64>,

    /// [DilutionOfPrecision]
    pub dop: DilutionOfPrecision,

    /// current [State]
    pub state: State,

    /// Null on first iter
    prev_epoch: Option<Epoch>,
}

impl Solver {
    /// Creates new [Solver].
    ///
    /// ## Input
    /// - cfg: [Config] preset
    /// - frame: [Frame]
    pub fn new(cfg: Config, frame: Frame) -> Self {
        let q_k = DMatrix::<f64>::zeros(U8::USIZE, U8::USIZE);
        let f_k = DMatrix::<f64>::identity(U8::USIZE, U8::USIZE);
        let g_k = DMatrix::<f64>::zeros(U8::USIZE, U8::USIZE);
        let w_k = DMatrix::<f64>::zeros(U8::USIZE, U8::USIZE);
        let p_k = DMatrix::<f64>::zeros(U8::USIZE, U8::USIZE);

        Self {
            f_k,
            q_k,
            g_k,
            w_k,
            p_k,
            cfg,
            frame,
            prev_epoch: None,
            state: Default::default(),
            sv: Vec::with_capacity(8),
            kalman: Kalman::new(U8::USIZE),
            y_k_vec: Vec::with_capacity(8),
            w_k_vec: Vec::with_capacity(8),
            indexes: Vec::with_capacity(8),
            fixed_amb: Default::default(),
            lambda_x: DMatrix::zeros(1, U4::USIZE),
            lambda_q: DMatrix::zeros(U4::USIZE, U4::USIZE),
            x_k: DVector::zeros(U8::USIZE),
            dop: DilutionOfPrecision::default(),
        }
    }

    /// Reset filter.
    pub fn reset(&mut self) {
        self.clear();
        self.fixed_amb.clear();

        self.kalman.reset();
        self.prev_epoch = None;
        self.dop = DilutionOfPrecision::default();
    }

    /// Iterates mutable [Navigation] filter.
    ///
    /// ## Input
    /// - epoch: sampling [Epoch]
    /// - params: [UserParameters]
    /// - candidates: proposed [Candidate]s
    /// - size: number of proposed [Cadndidate]s
    /// - rtk_base: custom [RTK] implementation
    /// - pivot_position_ecef_m: pivot position in meters (ECEF)
    /// - double_differences: double [Differences]
    pub fn run<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        initial_state: &State,
        candidates: &[Candidate],
        size: usize,
        rtk_base: &RTK,
        pivot_position_ecef_m: (f64, f64, f64),
        double_differences: &Differences,
    ) -> Result<(), Error> {
        self.clear();

        let initial_state = initial_state.clone();

        let mut ndf = U4::USIZE + double_differences.ndf();
        ndf -= 1; // TODO: only in RTK

        self.state.resize_mut(ndf);
        self.kalman.resize_mut(ndf);

        self.f_k.resize_mut(ndf, ndf, 0.0);
        self.q_k.resize_mut(ndf, ndf, 0.0);
        self.x_k.resize_vertically_mut(ndf, 0.0);

        if !self.kalman.initialized {
            self.kf_initialization(
                epoch,
                &initial_state,
                candidates,
                params,
                size,
                rtk_base,
                pivot_position_ecef_m,
                double_differences,
            )?;
        } else {
            self.kf_run(
                epoch,
                candidates,
                params,
                size,
                rtk_base,
                pivot_position_ecef_m,
                double_differences,
            )?;
        }

        self.prev_epoch = Some(epoch);

        Ok(())
    }

    /// Resets [Self] before a new run
    fn clear(&mut self) {
        self.sv.clear();
        self.indexes.clear();
        self.y_k_vec.clear();
        self.w_k_vec.clear();
    }

    /// Filter first iteration.
    ///
    /// ## Input
    /// - epoch: sampling [Epoch]
    /// - past_state: past [State]
    /// - candidates: proposed [Candidate]s
    /// - params: [UserParameters]
    /// - size: number of proposed [Cadndidate]s
    /// - rtk_base: custom [RTK] implementation
    /// - pivot_position_ecef_m: pivot position in meters (ECEF)
    /// - double_differences: double [Differences]
    pub fn kf_initialization<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        state: &State,
        candidates: &[Candidate],
        params: UserParameters,
        size: usize,
        rtk_base: &RTK,
        pivot_position_ecef_m: (f64, f64, f64),
        double_differences: &Differences,
    ) -> Result<(), Error> {
        const NB_ITER: usize = 10;

        let mut pending = state.clone();
        let mut dop = DilutionOfPrecision::default();

        let (base_x0, base_y0, base_z0) = rtk_base.reference_position_ecef_m(epoch);

        // measurements
        for i in 0..size {
            let mut contrib = SVContribution::default();

            contrib.sv = candidates[i].sv;

            match candidates[i].rtk_vector_contribution(
                epoch,
                true,
                &self.cfg,
                double_differences,
                &mut contrib,
            ) {
                Ok(vec) => {
                    self.y_k_vec.push(vec.row_1);
                    self.y_k_vec.push(vec.row_2);
                    self.w_k_vec.push(1.0); // TODO: improve model
                    self.w_k_vec.push(1.0); // TODO: improve model
                    self.indexes.push(i);
                    self.sv.push(contrib);
                },
                Err(e) => {
                    error!(
                        "{}({}) - ppp measurement error: {}",
                        epoch, candidates[i].sv, e
                    );
                },
            }
        }

        let y_len = self.y_k_vec.len();

        if y_len < U8::USIZE {
            return Err(Error::MatrixMinimalDimension);
        }

        // run
        for ith in 0..NB_ITER {
            let y_len = self.y_k_vec.len();

            // Build W
            self.w_k.resize_mut(y_len, y_len, 0.0);

            for i in 0..y_len {
                self.w_k[(i, i)] = 1.0 / self.w_k_vec[i];
            }

            let y_k = DVector::from_row_slice(&self.y_k_vec); // TODO malloc

            debug!("(ppp i={ith}) Y: {y_k}");

            // Build G
            let mut ndf = U4::USIZE;
            ndf -= 1; // TODO: only in RTK

            let lambda_ndf = self.indexes.len();
            ndf += lambda_ndf;
            self.state.resize_ambiguities_mut(lambda_ndf);

            self.g_k.resize_mut(y_len, ndf, 0.0);

            // Build G
            for (i, index) in self.indexes.iter().enumerate() {
                let position_m = pending.to_position_ecef_m();

                let (dx, dy, dz) =
                    candidates[*index].rtk_matrix_contribution(position_m, pivot_position_ecef_m);

                self.g_k[(2 * i, 0)] = dx;
                self.g_k[(2 * i, 1)] = dy;
                self.g_k[(2 * i, 2)] = dz;

                self.g_k[(2 * i + 1, 0)] = dx;
                self.g_k[(2 * i + 1, 1)] = dy;
                self.g_k[(2 * i + 1, 2)] = dz;
                self.g_k[(2 * i + 1, 3 + i)] = 1.0;
            }

            debug!("(ppp i={}) G: {}", ith, self.g_k);

            // run
            let gt = self.g_k.transpose();
            let gt_g = gt.clone() * self.g_k.clone();
            let gt_w = gt.clone() * self.w_k.clone();
            let gt_w_g = gt_w * self.g_k.clone();
            let gt_w_g_inv = gt_w_g.try_inverse().ok_or(Error::MatrixInversion)?;
            let gt_w_g_inv_gt = gt_w_g_inv.clone() * gt.clone();
            let gt_w_g_inv_gt_w = gt_w_g_inv_gt * self.w_k.clone();

            self.x_k = gt_w_g_inv_gt_w * y_k.clone();
            self.p_k = gt_w_g_inv.clone();

            let position_ecef_m = pending.to_position_ecef_m();

            let (baseline_dx, baseline_dy, baseline_dz) = (
                position_ecef_m[0] - base_x0,
                position_ecef_m[1] - base_y0,
                position_ecef_m[2] - base_z0,
            );

            self.x_k[0] -= baseline_dx;
            self.x_k[1] -= baseline_dy;
            self.x_k[2] -= baseline_dz;

            debug!("(ppp i={}) dx={}", ith, self.x_k);

            let (dx, dy, dz) = (self.x_k[0], self.x_k[1], self.x_k[2]);

            pending
                .spatial_correction_mut(self.frame, (dx, dy, dz))
                .map_err(|e| {
                    error!("{epoch} - state update failed with physical error: {e}");
                    Error::StateUpdate
                })?;

            pending.resize_ambiguities_mut(self.x_k.nrows() - 3); // TODO: only in RTK

            // ambiguities update
            for i in 3..self.x_k.nrows() {
                pending.x_amb[i - 3] += self.x_k[i];
            }

            let gt_g_inv = gt_g.try_inverse().ok_or(Error::MatrixInversion)?;

            // update latest DoP
            dop = DilutionOfPrecision::new(&pending, gt_g_inv);

            debug!("(ppp i={ith}) {epoch} - pending state {pending}");

            // models update
            self.y_k_vec.clear();
            self.w_k_vec.clear();

            self.indexes.retain(|i| {
                let mut unused = SVContribution::default();

                match candidates[*i].rtk_vector_contribution(
                    epoch,
                    true,
                    &self.cfg,
                    double_differences,
                    &mut unused,
                ) {
                    Ok(vec) => {
                        self.y_k_vec.push(vec.row_1);
                        self.y_k_vec.push(vec.row_2);
                        self.w_k_vec.push(1.0); // TODO improve model
                        self.w_k_vec.push(1.0); // TODO improve model
                        true
                    },
                    Err(e) => {
                        error!(
                            "{}({}) - rtk measurement error: {}",
                            epoch, candidates[*i].sv, e
                        );

                        false
                    },
                }
            });
        }

        // validation
        self.state_validation(&dop)?;

        debug!("dx(ppp)={}", self.x_k);

        let initial_estimate = KfEstimate::new(&self.x_k, &self.p_k);

        let mut ndf = U4::USIZE;
        ndf -= 1; // TODO: only in RTK

        let lambda_ndf = self.indexes.len();

        // TODO malloc
        let mut q_mat = DMatrix::identity(ndf + lambda_ndf, ndf + lambda_ndf);

        params.q_matrix(&mut q_mat, Duration::ZERO, ndf + lambda_ndf);

        for i in ndf..ndf + lambda_ndf {
            // TODO: code bias weight model
            // TODO: phase bias weight model
            q_mat[(i, i)] = 1.0;
        }

        debug!("ndf(ppp)={ndf} Q(ppp)={q_mat}");

        // build F
        self.f_k.resize_mut(ndf + lambda_ndf, ndf + lambda_ndf, 0.0);

        for i in 0..ndf + lambda_ndf {
            self.f_k[(i, i)] = 1.0;
        }

        self.kalman.initialize(&self.f_k, q_mat, initial_estimate);

        self.state = pending;

        debug!("{} - new PPP state {}", epoch, self.state);
        debug!("{} - gdop={} tdop={}", epoch, self.dop.gdop, self.dop.tdop);

        self.lambda_x.resize_mut(lambda_ndf, 1, 0.0);
        self.lambda_q.resize_mut(lambda_ndf, lambda_ndf, 0.0);

        let mut offset = U4::USIZE;
        offset -= 1; // TODO: only in RTK

        for i in 0..lambda_ndf {
            self.lambda_x[i] = self.x_k[offset + i];

            for j in 1..lambda_ndf {
                self.lambda_q[(i, j)] = self.p_k[(offset + i, j + offset)];
            }
        }

        match LambdaAR::run(lambda_ndf, lambda_ndf, &self.lambda_x, &self.lambda_q) {
            Ok((_f_mat, _s)) => {
                // TODO fix confirmation
                // TODO manque l'info de SV
                for (i, index) in self.indexes.iter().enumerate() {
                    let sv = candidates[*index].sv;
                    // TODO this is round() simplification (REMOVE)
                    self.fixed_amb
                        .insert(sv, self.x_k[offset + i].round() as i64);
                    // self.fixed_amb.insert(sv, f_mat[(i, 0)].round() as i64);
                }
            },
            Err(e) => {
                error!("{epoch} - lambda search failed with {e}");
                return Err(e);
            },
        }

        Ok(())
    }

    /// KF run
    ///
    /// ## Input
    /// - epoch: sampling [Epoch]
    /// - candidates: proposed [Candidate]s
    /// - params: [UserParameters]
    /// - size: number of proposed [Cadndidate]s
    /// - rtk_base: custom [RTK] implementation
    /// - pivot_position_ecef_m: pivot position in meters (ECEF)
    /// - double_differences: double [Differences]
    pub fn kf_run<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        candidates: &[Candidate],
        params: UserParameters,
        size: usize,
        rtk_base: &RTK,
        pivot_position_ecef_m: (f64, f64, f64),
        double_differences: &Differences,
    ) -> Result<(), Error> {
        let mut pending = self.state.clone();

        let (base_x0, base_y0, base_z0) = rtk_base.reference_position_ecef_m(epoch);

        // measurements
        for i in 0..size {
            let mut contrib = SVContribution::default();

            contrib.sv = candidates[i].sv;

            match candidates[i].rtk_vector_contribution(
                epoch,
                true,
                &self.cfg,
                double_differences,
                &mut contrib,
            ) {
                Ok(vec) => {
                    self.y_k_vec.push(vec.row_1);
                    self.y_k_vec.push(vec.row_2);
                    self.w_k_vec.push(1.0); // TODO
                    self.sv.push(contrib);
                    self.indexes.push(i);
                },
                Err(e) => {
                    error!(
                        "{}({}) - rtk measurement error: {}",
                        epoch, candidates[i].sv, e
                    );
                },
            }
        }

        let y_len = self.y_k_vec.len();

        if y_len < U8::USIZE {
            return Err(Error::MatrixMinimalDimension);
        }

        let y_k = DVector::from_row_slice(&self.y_k_vec); // TODO malloc
        debug!("Y(ppp): {y_k}");

        self.w_k.resize_mut(y_len, y_len, 0.0);

        let mut ndf = U4::USIZE;
        ndf -= 1; // TODO: only in RTK

        let lambda_ndf = self.indexes.len();

        // form W
        self.w_k.resize_mut(y_len, y_len, 0.0);

        for i in 0..y_len {
            self.w_k[(i, i)] = 1.0; // TODO: improve model
        }

        // form G
        self.g_k.resize_mut(y_len, ndf + lambda_ndf, 0.0);

        for (i, index) in self.indexes.iter().enumerate() {
            let position_m = pending.to_position_ecef_m();

            let (dx, dy, dz) =
                candidates[*index].rtk_matrix_contribution(position_m, pivot_position_ecef_m);

            self.g_k[(2 * i, 0)] = dx;
            self.g_k[(2 * i, 1)] = dy;
            self.g_k[(2 * i, 2)] = dz;

            self.g_k[(2 * i + 1, 0)] = dx;
            self.g_k[(2 * i + 1, 1)] = dy;
            self.g_k[(2 * i + 1, 2)] = dz;
        }

        debug!("G(ppp): {} W(ppp): {}", self.g_k, self.w_k);

        // form Y
        let y = DVector::<f64>::from_row_slice(self.y_k_vec.as_slice());

        debug!("Y(ppp): {y}");

        // Build F
        self.f_k = DMatrix::identity(ndf, ndf);

        // Build Q
        let mut q_mat = DMatrix::identity(ndf, ndf);

        params.q_matrix(&mut q_mat, Duration::ZERO, ndf);

        for i in ndf..ndf + lambda_ndf {
            // TODO: code bias weight model
            // TODO: phase bias weight model
            q_mat[(i, i)] = 1.0;
        }

        debug!("ndf(ppp)={} F(ppp)={} Q(ppp)={}", ndf, self.f_k, q_mat);

        let estimate = self
            .kalman
            .run(&self.f_k, &self.g_k, &self.w_k, &q_mat, &y)?;

        let ndf = estimate.x.nrows();

        let lambda_ndf = ndf - U6::USIZE; // TODO

        for i in 0..ndf {
            self.x_k[i] = estimate.x[i];
        }

        debug!("dx(ppp)={}", self.x_k);

        // if uses_rtk {
        let position_ecef_m = pending.to_position_ecef_m();

        let (baseline_dx, baseline_dy, baseline_dz) = (
            position_ecef_m[0] - base_x0,
            position_ecef_m[1] - base_y0,
            position_ecef_m[2] - base_z0,
        );

        self.x_k[0] -= baseline_dx;
        self.x_k[1] -= baseline_dy;
        self.x_k[2] -= baseline_dz;
        // }

        let (dx, dy, dz) = (self.x_k[0], self.x_k[1], self.x_k[2]);

        pending
            .spatial_correction_mut(self.frame, (dx, dy, dz))
            .map_err(|e| {
                error!("{epoch} - state update failed with physical error: {e}");
                Error::StateUpdate
            })?;

        let gt_g_inv = (self.g_k.transpose() * self.g_k.clone())
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        // DoP update
        let dop = DilutionOfPrecision::new(&pending, gt_g_inv);

        // validation
        self.state_validation(&dop)?;

        self.state = pending;
        self.dop = dop;

        debug!("{} - new PPP state {}", epoch, self.state);
        debug!("{} - gdop={} tdop={}", epoch, self.dop.gdop, self.dop.tdop);

        if self.cfg.method == Method::PPP {
            self.lambda_q.resize_mut(lambda_ndf, lambda_ndf, 0.0);
            self.lambda_x.resize_mut(lambda_ndf, 1, 0.0);

            for i in 0..lambda_ndf {
                for j in 0..lambda_ndf {
                    self.lambda_q[(i, j)] = self.p_k[(
                        i + Navigation::clock_index() + 1,
                        j + Navigation::clock_index() + 1,
                    )];
                }

                self.lambda_x[i] = self.x_k[i + Navigation::clock_index() + 1];
            }

            match LambdaAR::run(lambda_ndf, lambda_ndf, &self.lambda_x, &self.lambda_q) {
                Ok(_) => {
                    // TODO confirmation
                },
                Err(e) => {
                    error!("lambda search failed with {e}");
                    return Err(e);
                },
            }
        }

        Ok(())
    }

    /// Validate pending [State]
    fn state_validation(&self, dop: &DilutionOfPrecision) -> Result<(), Error> {
        if dop.gdop > self.cfg.solver.max_gdop {
            return Err(Error::MaxGdopExceeded);
        }

        Ok(())
    }
}
