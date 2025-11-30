use nalgebra::{DimName, Matrix6, OMatrix, U3, U6, Vector6};

use crate::{
    error::Error,
    navigation::{
        kalman::{Kalman, KfEstimate},
        state::State,
    },
    prelude::Duration,
};

#[derive(Clone)]
pub struct PostfitKf {
    /// F [OMatrix]
    f_mat: OMatrix<f64, U6, U6>,

    /// G [OMatrix]
    g_mat: OMatrix<f64, U6, U6>,

    /// W [OMatrix]
    w_mat: OMatrix<f64, U6, U6>,

    /// Q [OMatrix]
    q_mat: OMatrix<f64, U6, U6>,

    /// [Kalman] filter
    kalman: Kalman,
}

impl PostfitKf {
    /// Builds new [PostfitKf] from initial [State]
    pub fn new(
        state: &State,
        state_pos_std_dev_m: f64,
        state_vel_std_dev_m_s: f64,
        meas_pos_std_dev_m: f64,
        meas_vel_std_dev_m_s: f64,
    ) -> Self {
        let x_0 = state.to_position_velocity_ecef_m();

        let q_diag = Vector6::new(
            state_pos_std_dev_m.powi(2),
            state_pos_std_dev_m.powi(2),
            state_pos_std_dev_m.powi(2),
            state_vel_std_dev_m_s.powi(2),
            state_vel_std_dev_m_s.powi(2),
            state_vel_std_dev_m_s.powi(2),
        );

        let q_mat = Matrix6::from_diagonal(&q_diag);
        let p_0 = Matrix6::from_diagonal(&q_diag);
        let f_mat = Matrix6::identity();

        let mut kalman = Kalman::new(U6::DIM);

        let initial_estimate = KfEstimate::from_static(x_0, p_0);

        kalman.initialize_from_static(f_mat, q_mat, initial_estimate);

        let mut w_mat = Matrix6::identity();

        for i in 0..U3::DIM {
            w_mat[(i, i)] = meas_pos_std_dev_m;
            w_mat[(i + U3::DIM, i + U3::DIM)] = meas_vel_std_dev_m_s;
        }

        Self {
            q_mat,
            kalman,
            f_mat: Matrix6::identity(),
            g_mat: Matrix6::identity(),
            w_mat: Matrix6::identity(),
        }
    }

    /// Run [PostfitKf] filter.
    ///
    /// ## Input
    /// - state: new [State]
    /// - sampling_interval: [Duration]
    ///
    /// ## Output
    /// - estimate: [KfEstimate]
    pub fn run(&mut self, state: &State, sampling_interval: Duration) -> Result<KfEstimate, Error> {
        let dt_s = sampling_interval.to_seconds();

        self.f_mat[(0, 3)] = dt_s;
        self.f_mat[(1, 4)] = dt_s;
        self.f_mat[(2, 5)] = dt_s;

        let y_vec = &state.to_position_velocity_ecef_m();

        self.kalman
            .run_static::<U6>(&self.f_mat, &self.g_mat, &self.w_mat, &self.q_mat, y_vec)
    }

    /// Reset this [PostfitKf]
    pub fn reset(&mut self) {
        self.kalman.reset();
    }
}
