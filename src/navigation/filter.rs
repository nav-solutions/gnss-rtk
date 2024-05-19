use nalgebra::{base::dimension::U8, OMatrix, OVector, Vector3};

#[cfg(feature = "serde")]
use serde::Deserialize;

use super::{Input, Output};
use crate::prelude::{Epoch, Error};

/// Navigation Filter.
#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum Filter {
    /// None: solver filter completely bypassed. Lighter calculations, no iterative behavior.
    None,
    #[default]
    /// LSQ Filter. Heavy computation.
    LSQ,
    /// Kalman Filter. Heavy+ computations. Compared to LSQ, the Kalman filter
    /// converges faster and has the ability to "improve" models
    Kalman,
}

#[derive(Debug, Clone, Default)]
struct LSQState {
    pub p: OMatrix<f64, U8, U8>,
    pub x: OVector<f64, U8>,
}

#[derive(Debug, Clone, Default)]
struct KFState {
    pub q: OMatrix<f64, U8, U8>,
    pub p: OMatrix<f64, U8, U8>,
    pub x: OVector<f64, U8>,
    pub phi: OMatrix<f64, U8, U8>,
}

#[derive(Debug, Clone)]
pub enum FilterState {
    Lsq(LSQState),
    Kf(KFState),
}

impl Default for FilterState {
    fn default() -> Self {
        Self::Lsq(Default::default())
    }
}

impl FilterState {
    fn lsq(state: LSQState) -> Self {
        Self::Lsq(state)
    }
    //fn as_lsq(&self) -> Option<&LSQState> {
    //    match self {
    //        Self::LSQ(state) => Some(state),
    //        _ => None,
    //    }
    //}
    pub fn ambiguities(&self) -> Vec<f64> {
        let x = match self {
            Self::Lsq(state) => state.x,
            Self::Kf(state) => state.x,
        };
        let mut r = Vec::<f64>::new();
        for i in 4..x.len() {
            r.push(x[i]);
        }
        r
    }
    fn kf(state: KFState) -> Self {
        Self::Kf(state)
    }
    //fn as_kf(&self) -> Option<&KFState> {
    //    match self {
    //        Self::KF(state) => Some(state),
    //        _ => None,
    //    }
    //}
    pub(crate) fn estimate(&self) -> OVector<f64, U8> {
        match self {
            Self::Lsq(state) => state.x,
            Self::Kf(state) => state.x,
        }
    }
}

impl Filter {
    fn lsq_resolve(input: &Input, p_state: Option<FilterState>) -> Result<Output, Error> {
        match p_state {
            Some(FilterState::Lsq(p_state)) => {
                let p_1 = p_state.p.try_inverse().ok_or(Error::MatrixInversionError)?;

                let g_prime = input.g.clone().transpose();
                let q = (g_prime * input.g)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = g_prime * input.w * input.g;
                let p = (p_1 + p).try_inverse().ok_or(Error::MatrixInversionError)?;

                let x = p * (p_1 * p_state.x + (g_prime * input.w * input.y));

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::lsq(LSQState { p, x }),
                })
            },
            _ => {
                let g_prime = input.g.clone().transpose();

                let q = (g_prime * input.g)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = (g_prime * input.w * input.g)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let x = p * (g_prime * input.w * input.y);
                if x[3].is_nan() {
                    return Err(Error::TimeIsNan);
                }

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::lsq(LSQState { p, x }),
                })
            },
        }
    }
    fn kf_resolve(input: &Input, p_state: Option<FilterState>) -> Result<Output, Error> {
        match p_state {
            Some(FilterState::Kf(p_state)) => {
                let x_bn = p_state.phi * p_state.x;
                let p_bn = p_state.phi * p_state.p * p_state.phi.transpose() + p_state.q;

                let p_bn_inv = p_bn.try_inverse().ok_or(Error::MatrixInversionError)?;
                let p_n = (input.g.transpose() * input.w * input.g + p_bn_inv)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let w_g = input.g.transpose() * input.w * input.y;
                let w_gy_pbn = w_g + (p_bn_inv * x_bn);
                let x_n = p_n * w_gy_pbn;

                let q_n = input.g.transpose() * input.g;
                let phi_diag = OVector::<f64, U8>::from([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
                let q_diag = OVector::<f64, U8>::from([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);

                Ok(Output {
                    gdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)] + q_n[(3, 3)]).sqrt(),
                    pdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)]).sqrt(),
                    tdop: q_n[(4, 3)].sqrt(),
                    q: q_n,
                    state: FilterState::kf(KFState {
                        p: p_n,
                        x: x_n,
                        q: OMatrix::<f64, U8, U8>::from_diagonal(&q_diag),
                        phi: OMatrix::<f64, U8, U8>::from_diagonal(&phi_diag),
                    }),
                })
            },
            _ => {
                let g_prime = input.g.clone().transpose();
                let q = (g_prime * input.g)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = (g_prime * input.w * input.g)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let x = p * (g_prime * input.w * input.y);
                if x[3].is_nan() {
                    return Err(Error::TimeIsNan);
                }

                let phi_diag = OVector::<f64, U8>::from([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
                let q_diag = OVector::<f64, U8>::from([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::kf(KFState {
                        p,
                        x,
                        q: OMatrix::<f64, U8, U8>::from_diagonal(&q_diag),
                        phi: OMatrix::<f64, U8, U8>::from_diagonal(&phi_diag),
                    }),
                })
            },
        }
    }
    pub fn resolve(&self, input: &Input, p_state: Option<FilterState>) -> Result<Output, Error> {
        match self {
            Filter::None => Self::lsq_resolve(input, None),
            Filter::LSQ => Self::lsq_resolve(input, p_state),
            Filter::Kalman => Self::kf_resolve(input, p_state),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Copy, Default)]
pub(crate) struct State3D {
    pub t: Epoch,
    pub inner: Vector3<f64>,
}

impl std::fmt::LowerExp for State3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let (x, y, z) = (self.inner[0], self.inner[1], self.inner[2]);
        write!(f, "({:.6E},{:.6E},{:.6E})", x, y, z)
    }
}

impl std::fmt::Display for State3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let (x, y, z) = (self.inner[0], self.inner[1], self.inner[2]);
        write!(f, "({:.6E},{:.6E},{:.6E})", x, y, z)
    }
}

// impl NyxState for State3D {
//     type Size = U3;
//     type VecLength = U3;
//     fn as_vector(&self) -> OVector<f64, U3>, NyxError {
//         self.inner.into()
//     }
//     fn unset_stm(&mut self) {}
//     fn set(&mut self, t: Epoch, vector: &OVector<f64, U3>) -> () {
//         self.t = t;
//         self.inner = vector.clone();
//     }
//     fn epoch(&self) -> Epoch {
//         self.t
//     }
//     fn set_epoch(&mut self, t: Epoch) {
//         self.t = t;
//     }
// }
