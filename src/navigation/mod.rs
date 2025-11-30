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

// use std::collections::HashMap;

use nalgebra::{DMatrix, DVector, DimName, U3, U4};

use crate::{
    candidate::differences::Differences,
    navigation::{
        apriori::Apriori,
        dop::DilutionOfPrecision,
        kalman::{Kalman, KfEstimate},
        postfit::PostfitKf,
        ppp_ar::Solver,
        state::State,
        sv::SVContribution,
    },
    prelude::{
        Candidate,
        Config,
        Duration,
        Epoch,
        Error,
        Frame,
        // Method,
        UserParameters,
        // SV,
    },
    rtk::RTKBase,
};

/// [Navigation] Solver
pub(crate) struct Navigation {
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

    /// G
    g_k: DMatrix<f64>,

    /// F
    f_k: DMatrix<f64>,

    /// Q
    q_k: DMatrix<f64>,

    /// X
    x_k: DVector<f64>,

    /// P
    p_k: DMatrix<f64>,

    /// Prefit
    prefit: Option<Solver>,

    /// [Kalman]
    kalman: Kalman,

    /// Postfit (Kalman)
    postfit: Option<PostfitKf>,

    /// contribution allocation
    indexes: Vec<usize>,

    /// [SVContribution]s
    pub sv: Vec<SVContribution>,

    /// Current [State]
    pub state: State,

    /// [DilutionOfPrecision]
    pub dop: DilutionOfPrecision,

    /// Null on first iter
    prev_epoch: Option<Epoch>,
}

impl Navigation {
    /// Creates new [Navigation] solver.
    ///
    /// ## Input
    /// - cfg: [Config] preset
    /// - apriori: [Apriori] input
    /// - candidates: selected [Candidate]s
    /// - size: number of proposal
    /// ## Returns
    /// - [Navigation], [Error]
    pub fn new(cfg: &Config, frame: Frame) -> Self {
        let q_k = DMatrix::<f64>::zeros(U4::DIM, U4::DIM);
        let mut f_k = DMatrix::<f64>::zeros(U4::DIM, U4::DIM);

        for i in 0..=U3::DIM {
            f_k[(i, i)] = 1.0;
        }

        Self {
            f_k,
            q_k,
            frame,
            postfit: None,
            prev_epoch: None,
            cfg: cfg.clone(),
            prefit: None,
            state: Default::default(),
            sv: Vec::with_capacity(8),
            kalman: Kalman::new(U4::DIM),
            y_k_vec: Vec::with_capacity(8),
            w_k_vec: Vec::with_capacity(8),
            indexes: Vec::with_capacity(8),
            x_k: DVector::zeros(U4::DIM),
            dop: DilutionOfPrecision::default(),
            g_k: DMatrix::<f64>::zeros(U4::DIM, U4::DIM),
            w_k: DMatrix::<f64>::zeros(U4::DIM, U4::DIM),
            p_k: DMatrix::<f64>::zeros(U4::DIM, U4::DIM),
        }
    }

    /// Reset [Navigation] filter.
    ///
    /// To this day, it has not been demonstrated that this
    /// operation is needed by any real-world application, even
    /// real-time processing. But anyways.. this exists.
    ///
    /// After reset, the [Navigation] solver is in the same state as deployment,
    /// ready to consume a first measurement.
    pub fn reset(&mut self) {
        self.clear();
        self.kalman.reset();

        if let Some(prefit) = &mut self.prefit {
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
    /// - user: [User] profile
    /// - past_state: past [State]
    /// - candidates: proposed [Candidate]s
    /// - size: number of proposed [Cadndidate]s
    /// - uses_rtk: true when RTK mode nav is being used
    /// - double_differences: possible double [Differences]
    pub fn solve<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        initial_state: &State,
        candidates: &[Candidate],
        size: usize,
        uses_rtk: bool,
        rtk_base: &RTK,
        pivot_position_ecef_m: &Option<(f64, f64, f64)>,
        double_differences: &Option<Differences>,
    ) -> Result<(), Error> {
        self.clear();

        let initial_state = initial_state.clone();

        let mut ndf = U4::DIM;

        if uses_rtk {
            ndf -= 1;
        }

        self.state.resize_mut(ndf);
        self.kalman.resize_mut(ndf);

        self.f_k.resize_mut(ndf, ndf, 0.0);
        self.q_k.resize_mut(ndf, ndf, 0.0);
        self.x_k.resize_vertically_mut(ndf, 0.0);

        for i in 0..ndf {
            self.f_k[(i, i)] = 1.0;
        }

        let dt = if let Some(past_epoch) = self.prev_epoch {
            epoch - past_epoch
        } else {
            Duration::ZERO
        };

        // if self.cfg.method == Method::PPP {
        //     if self.prefit.is_none() {
        //         self.prefit = Some(Solver::new(self.cfg.clone(), self.frame));
        //     }
        // }

        params.q_matrix(&mut self.q_k, dt, ndf);

        if let Some(prefit) = &mut self.prefit {
            let double_diff = double_differences
                .as_ref()
                .expect("internal error: missing RTK+PPP prefit");

            let pivot_position_ecef_m = pivot_position_ecef_m.unwrap_or_else(|| {
                panic!("internal error: undefined pivot satellite position");
            });

            match prefit.run(
                epoch,
                params,
                &initial_state,
                candidates,
                size,
                rtk_base,
                pivot_position_ecef_m,
                double_diff,
            ) {
                Ok(_) => {},
                Err(e) => {
                    error!("{epoch} - ppp prefit failed with {e}");
                    return Err(Error::FloatAmbiguitiesSolving);
                },
            }
        }

        // let fixed_ambiguities = if self.cfg.method == Method::PPP {
        //     // let prefit = self
        //     //     .prefit
        //     //     .as_ref()
        //     //     .expect("internal error: missing prefit solver");

        //     // Some(prefit.fixed_amb.clone())
        //     None
        // } else {
        //     None
        // };

        if !self.kalman.initialized {
            self.kf_initialization(
                epoch,
                &initial_state,
                candidates,
                size,
                uses_rtk,
                rtk_base,
                pivot_position_ecef_m,
                double_differences,
            )?;
        } else {
            self.kf_run(
                epoch,
                candidates,
                size,
                uses_rtk,
                rtk_base,
                pivot_position_ecef_m,
                double_differences,
            )?;
        }

        if self.cfg.solver.postfit_denoising > 0.0 {
            if let Some(postfit) = &mut self.postfit {
                let prev_epoch = self
                    .prev_epoch
                    .expect("internal error: undetermined past epoch");

                let dt = epoch - prev_epoch;
                let dx = postfit.run(&self.state, dt)?;

                // update state
                self.state
                    .postfit_update_mut(self.frame, &dx.x)
                    .map_err(|e| {
                        error!("{epoch} - postfit state update failed with physical error: {e}");
                        Error::StateUpdate
                    })?;
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
    /// - user: [User] profile
    /// - past_state: past [State]
    /// - candidates: proposed [Candidate]s
    /// - size: number of proposed [Cadndidate]s
    /// - uses_rtk: true when RTK mode nav is being used
    /// - double_differences: possible double [Differences]
    /// - fixed_ambiguities: possible fixed ambiguities
    pub fn kf_initialization<RTK: RTKBase>(
        &mut self,
        t: Epoch,
        state: &State,
        candidates: &[Candidate],
        size: usize,
        uses_rtk: bool,
        rtk_base: &RTK,
        pivot_position_ecef_m: &Option<(f64, f64, f64)>,
        double_differences: &Option<Differences>,
    ) -> Result<(), Error> {
        const NB_ITER: usize = 10;

        let mut pending = state.clone();
        let mut dop = DilutionOfPrecision::default();

        let (base_x0, base_y0, base_z0) = rtk_base.reference_position_ecef_m(t);

        // measurements
        for i in 0..size {
            let mut contrib = SVContribution::default();

            contrib.sv = candidates[i].sv;

            let position_m = pending.to_position_ecef_m();

            let vec = if uses_rtk {
                let double_differences = double_differences
                    .as_ref()
                    .expect("internal error: invalid rtk measurement/post fit");

                // let fixed_amb = if self.cfg.method == Method::PPP {
                //     let fixed_ambiguities = fixed_ambiguities
                //         .as_ref()
                //         .expect("internal error: missing PPP prefit");

                //     fixed_ambiguities.get(&candidates[i].sv)
                // } else {
                //     None
                // };

                match candidates[i].rtk_vector_contribution(
                    t,
                    false,
                    &self.cfg,
                    double_differences,
                    &mut contrib,
                ) {
                    Ok(vec) => Some(vec),
                    Err(e) => {
                        error!("{}({}) - rtk measurement error: {}", t, candidates[i].sv, e);
                        None
                    },
                }
            } else {
                match candidates[i].ppp_vector_contribution(
                    &self.cfg,
                    false,
                    position_m,
                    &mut contrib,
                ) {
                    Ok(vec) => Some(vec),
                    Err(e) => {
                        error!("{}({}) - ppp measurement error: {}", t, candidates[i].sv, e);
                        None
                    },
                }
            };

            if vec.is_none() {
                continue;
            }

            let vec = vec.unwrap();

            self.y_k_vec.push(vec.row_1);
            self.w_k_vec.push(1.0); // TODO improve model
            self.indexes.push(i);
            self.sv.push(contrib);
        }

        let y_len = self.indexes.len();

        if y_len < U4::DIM {
            return Err(Error::MatrixMinimalDimension);
        }

        let mut ndf = U4::DIM;

        if uses_rtk {
            ndf -= 1;
        }

        self.g_k.resize_mut(y_len, ndf, 0.0);

        // run
        for ith in 0..NB_ITER {
            let y_len = self.y_k_vec.len();

            // Form W
            self.w_k.resize_mut(y_len, y_len, 0.0);

            for i in 0..y_len {
                self.w_k[(i, i)] = 1.0 / self.w_k_vec[i];
            }

            let y_k = DVector::from_row_slice(&self.y_k_vec); // TODO malloc
            debug!("(i={ith}) Y: {y_k}");

            // Form G
            for (i, index) in self.indexes.iter().enumerate() {
                let position_m = pending.to_position_ecef_m();

                let (dx, dy, dz) = if uses_rtk {
                    let pivot_position_ecef_m = pivot_position_ecef_m.unwrap_or_else(|| {
                        panic!("internal error: undefined pivot satellite position");
                    });

                    candidates[*index].rtk_matrix_contribution(position_m, pivot_position_ecef_m)
                } else {
                    candidates[*index].ppp_matrix_contribution(&self.cfg, position_m)
                };

                self.g_k[(i, 0)] = dx;
                self.g_k[(i, 1)] = dy;
                self.g_k[(i, 2)] = dz;

                if !uses_rtk {
                    self.g_k[(i, Self::clock_index())] = 1.0;
                }
            }

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

            if uses_rtk {
                let position_ecef_m = pending.to_position_ecef_m();

                let (baseline_dx, baseline_dy, baseline_dz) = (
                    position_ecef_m[0] - base_x0,
                    position_ecef_m[1] - base_y0,
                    position_ecef_m[2] - base_z0,
                );

                self.x_k[0] -= baseline_dx;
                self.x_k[1] -= baseline_dy;
                self.x_k[2] -= baseline_dz;
            }

            debug!("(i={}) dx={}", ith, self.x_k);

            let (dx, dy, dz) = (self.x_k[0], self.x_k[1], self.x_k[2]);

            pending
                .spatial_correction_mut(self.frame, (dx, dy, dz))
                .map_err(|e| {
                    error!("{t} - state update failed with physical error: {e}");
                    Error::StateUpdate
                })?;

            if !uses_rtk {
                pending.temporal_correction_mut(self.x_k[3]);
            }

            let gt_g_inv = gt_g.try_inverse().ok_or(Error::MatrixInversion)?;

            // update latest DoP
            dop = DilutionOfPrecision::new(&pending, gt_g_inv);

            debug!("(i={ith}) {t} - pending state {pending}");

            // models update
            self.y_k_vec.clear();
            self.w_k_vec.clear();

            self.indexes.retain(|i| {
                let position_m = pending.to_position_ecef_m();

                let mut unused = SVContribution::default();

                let vec = if uses_rtk {
                    let double_differences = double_differences
                        .as_ref()
                        .expect("internal error: invalid rtk measurement/post fit");

                    // let fixed_amb = if self.cfg.method == Method::PPP {
                    //     let fixed_ambiguities = fixed_ambiguities
                    //         .as_ref()
                    //         .expect("internal error: missing PPP prefit");

                    //     fixed_ambiguities.get(&candidates[*i].sv)
                    // } else {
                    //     None
                    // };

                    match candidates[*i].rtk_vector_contribution(
                        t,
                        false,
                        &self.cfg,
                        double_differences,
                        &mut unused,
                    ) {
                        Ok(vec) => Some(vec),
                        Err(e) => {
                            error!(
                                "{}({}) - rtk measurement error: {}",
                                t, candidates[*i].sv, e
                            );
                            None
                        },
                    }
                } else {
                    match candidates[*i].ppp_vector_contribution(
                        &self.cfg,
                        false,
                        position_m,
                        &mut unused,
                    ) {
                        Ok(vec) => Some(vec),
                        Err(e) => {
                            error!(
                                "{}({}) - ppp measurement error: {}",
                                t, candidates[*i].sv, e
                            );
                            None
                        },
                    }
                };

                if let Some(vec) = vec {
                    self.y_k_vec.push(vec.row_1);
                    self.w_k_vec.push(1.0); // TODO improve model
                    true
                } else {
                    false
                }
            });
        }

        // validation
        self.state_validation(&dop)?;

        let initial_estimate = KfEstimate::new(&self.x_k, &self.p_k);

        self.kalman
            .initialize(&self.f_k, self.q_k.clone(), initial_estimate);

        self.state = pending;
        self.dop = dop;

        debug!("{} - new state {}", t, self.state);
        debug!("{} - gdop={} tdop={}", t, self.dop.gdop, self.dop.tdop);

        Ok(())
    }

    /// [Kalman] filter run
    ///
    /// ## Input
    /// - epoch: sampling [Epoch]
    /// - user: [User] profile
    /// - past_state: past [State]
    /// - candidates: proposed [Candidate]s
    /// - size: number of proposed [Cadndidate]s
    /// - uses_rtk: true when RTK mode nav is being used
    /// - double_differences: possible double [Differences]
    pub fn kf_run<RTK: RTKBase>(
        &mut self,
        t: Epoch,
        candidates: &[Candidate],
        size: usize,
        uses_rtk: bool,
        rtk_base: &RTK,
        pivot_position_ecef_m: &Option<(f64, f64, f64)>,
        double_differences: &Option<Differences>,
    ) -> Result<(), Error> {
        let mut pending = self.state.clone();

        let (base_x0, base_y0, base_z0) = rtk_base.reference_position_ecef_m(t);

        // measurement
        for i in 0..size {
            let mut contrib = SVContribution::default();

            contrib.sv = candidates[i].sv;

            let pos_m = pending.to_position_ecef_m();

            let vec = if uses_rtk {
                let double_differences = double_differences
                    .as_ref()
                    .expect("internal error: invalid rtk measurement/post fit");

                // let fixed_amb = if self.cfg.method == Method::PPP {
                //     // let fixed_ambiguities = fixed_ambiguities
                //     //     .as_ref()
                //     //     .expect("internal error: missing PPP prefit");

                //     // fixed_ambiguities.get(&candidates[i].sv)
                //     None
                // } else {
                //     None
                // };

                match candidates[i].rtk_vector_contribution(
                    t,
                    false,
                    &self.cfg,
                    double_differences,
                    &mut contrib,
                ) {
                    Ok(vec) => Some(vec),
                    Err(e) => {
                        error!("{}({}) - rtk measurement error: {}", t, candidates[i].sv, e);
                        None
                    },
                }
            } else {
                match candidates[i].ppp_vector_contribution(&self.cfg, false, pos_m, &mut contrib) {
                    Ok(vec) => Some(vec),
                    Err(e) => {
                        error!("{}({}) - ppp measurement error: {}", t, candidates[i].sv, e);
                        None
                    },
                }
            };

            if let Some(vec) = vec {
                self.y_k_vec.push(vec.row_1);
                self.w_k_vec.push(1.0); // TODO improve model
                self.indexes.push(i);
                self.sv.push(contrib);
            } else {
                error!("{}({}) - cannot contribute", t, candidates[i].sv);
            }
        }

        let y_len = self.y_k_vec.len();

        if y_len < U4::DIM {
            return Err(Error::MatrixMinimalDimension);
        }

        let y_k = DVector::from_row_slice(&self.y_k_vec); // TODO malloc
        debug!("Y: {y_k}");

        self.w_k.resize_mut(y_len, y_len, 0.0);

        let mut ndf = U4::DIM;

        if uses_rtk {
            ndf -= 1;
        }

        self.g_k.resize_mut(y_len, ndf, 0.0);

        for i in 0..y_len {
            self.w_k[(i, i)] = 1.0 / self.w_k_vec[i];
        }

        // Form G
        for (i, index) in self.indexes.iter().enumerate() {
            let position_m = pending.to_position_ecef_m();

            let (dx, dy, dz) = if uses_rtk {
                let pivot_position_ecef_m = pivot_position_ecef_m.unwrap_or_else(|| {
                    panic!("internal error: undefined pivot satellite position");
                });

                candidates[*index].rtk_matrix_contribution(position_m, pivot_position_ecef_m)
            } else {
                candidates[*index].ppp_matrix_contribution(&self.cfg, position_m)
            };

            self.g_k[(i, 0)] = dx;
            self.g_k[(i, 1)] = dy;
            self.g_k[(i, 2)] = dz;

            if !uses_rtk {
                self.g_k[(i, Self::clock_index())] = 1.0;
            }
        }

        let estimate = self
            .kalman
            .run(&self.f_k, &self.g_k, &self.w_k, &self.q_k, &y_k)?;

        debug!("state correction: dx={}", estimate.x);

        let ndf = estimate.x.nrows();

        for i in 0..ndf {
            self.x_k[i] = estimate.x[i];
        }

        if uses_rtk {
            let position_ecef_m = pending.to_position_ecef_m();

            let (baseline_dx, baseline_dy, baseline_dz) = (
                position_ecef_m[0] - base_x0,
                position_ecef_m[1] - base_y0,
                position_ecef_m[2] - base_z0,
            );

            self.x_k[0] -= baseline_dx;
            self.x_k[1] -= baseline_dy;
            self.x_k[2] -= baseline_dz;
        }

        let (dx, dy, dz) = (self.x_k[0], self.x_k[1], self.x_k[2]);

        pending
            .spatial_correction_mut(self.frame, (dx, dy, dz))
            .map_err(|e| {
                error!("{t} - state update failed with physical error: {e}");
                Error::StateUpdate
            })?;

        if !uses_rtk {
            pending.temporal_correction_mut(self.x_k[3]);
        }

        let gt_g_inv = (self.g_k.transpose() * self.g_k.clone())
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        // update
        let dop = DilutionOfPrecision::new(&pending, gt_g_inv);

        self.state_validation(&dop)?;

        self.dop = dop;
        self.state = pending.clone();

        debug!("{t} - new state {pending}");
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
    use crate::navigation::{DilutionOfPrecision, Navigation, State};
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
}
