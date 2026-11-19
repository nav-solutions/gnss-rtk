use crate::error::Error;

use nalgebra::{
    DMatrix, DVector, DefaultAllocator, DimName, OMatrix, OVector, allocator::Allocator,
};

#[derive(Clone)]
pub struct KfEstimate {
    /// P [DMatrix]
    pub p: DMatrix<f64>,

    /// x [DVector]
    pub x: DVector<f64>,
}

impl KfEstimate {
    /// Create a zero [KfEstimate]
    pub fn zero(size: usize) -> Self {
        let x = DVector::<f64>::zeros(size);
        let p = DMatrix::<f64>::zeros(size, size);
        Self { p, x }
    }

    pub fn resize_mut(&mut self, ndf: usize) {
        self.p.resize_mut(ndf, ndf, 100.0e-3);
        self.x.resize_vertically_mut(ndf, 0.0);
    }

    /// Reset this [KfEstimate]
    pub fn reset(&mut self) {
        for i in 0..self.x.nrows() {
            self.x[i] = 0.0;

            for j in 0..self.p.ncols() {
                self.p[(i, j)] = 0.0;
            }
        }
    }

    /// Initializes new [KfEstimate]
    pub fn new(x: &DVector<f64>, p: &DMatrix<f64>) -> Self {
        let x_rows = x.nrows();
        let (p_rows, p_cols) = (p.nrows(), p.ncols());

        assert_eq!(x_rows, p_cols, "P/x dimension issue!");
        assert_eq!(p_rows, p_cols, "P is not square!");

        Self {
            p: p.clone(),
            x: x.clone(),
        }
    }

    // /// Initializes [KfEstimate] from dyn/static matrices
    // pub fn from_dyn_static<D: DimName>(
    //     x: Matrix<f64, D, Dyn, VecStorage<f64, D, Dyn>>,
    //     p: OMatrix<f64, D, D>,
    // ) -> Self
    // where
    //     D: DimName,
    //     DefaultAllocator: nalgebra::allocator::Allocator<D>,
    //     DefaultAllocator: nalgebra::allocator::Allocator<D, D>,
    // {
    //     let (x_rows, x_cols) = (x.nrows(), x.ncols());
    //     let (p_rows, p_cols) = (p.nrows(), p.ncols());

    //     assert_eq!(x_cols, 1, "x is not a column vector!");
    //     assert_eq!(x_rows, p_cols, "P/x dimension issue!");
    //     assert_eq!(p_rows, p_cols, "P is not square!");

    //     Self {
    //         x: DMatrix::from_column_slice(D::DIM, U1::DIM, x.as_slice()),
    //         p: DMatrix::from_column_slice(D::DIM, D::DIM, p.as_slice()),
    //     }
    // }

    /// Initializes [KfEstimate] from static matrices
    pub fn from_static<D: DimName>(x: OVector<f64, D>, p: OMatrix<f64, D, D>) -> Self
    where
        D: DimName,
        DefaultAllocator: nalgebra::allocator::Allocator<D>,
        DefaultAllocator: nalgebra::allocator::Allocator<D, D>,
    {
        let (x_rows, _) = (x.nrows(), x.ncols());
        let (p_rows, p_cols) = (p.nrows(), p.ncols());

        assert_eq!(x_rows, p_cols, "P/x dimension issue!");
        assert_eq!(p_rows, p_cols, "P is not square!");

        Self {
            x: DVector::from_column_slice(x.as_slice()),
            p: DMatrix::from_column_slice(D::DIM, D::DIM, p.as_slice()),
        }
    }
}

#[derive(Clone)]
pub struct Kalman {
    /// True if this [Kalman] filter has been initialized
    pub initialized: bool,

    /// Prediction as [KfEstimate].
    pub predicted: KfEstimate,
}

impl Kalman {
    /// Create a new [Kalman] filter
    pub fn new(size: usize) -> Self {
        Self {
            initialized: false,
            predicted: KfEstimate::zero(size),
        }
    }

    pub fn resize_mut(&mut self, ndf: usize) {
        self.predicted.resize_mut(ndf);
    }

    /// Initialize this [Kalman] filter
    pub fn initialize_from_static<D: DimName>(
        &mut self,
        f_k: OMatrix<f64, D, D>,
        q_k: OMatrix<f64, D, D>,
        estimate: KfEstimate,
    ) where
        D: DimName,
        DefaultAllocator: Allocator<D>,
        DefaultAllocator: Allocator<D, D>,
    {
        // prediction
        let x_k = f_k.clone() * estimate.x;
        let f_k_t = f_k.transpose();
        let p_k = f_k * estimate.p * f_k_t + q_k;

        let x_k = DVector::<f64>::from_row_slice(x_k.as_slice());
        let p_k = DMatrix::<f64>::from_column_slice(D::DIM, D::DIM, p_k.as_slice());

        self.predicted = KfEstimate::new(&x_k, &p_k);
        self.initialized = true;
    }

    /// Initialize this [Kalman] filter
    pub fn initialize(&mut self, f_k: &DMatrix<f64>, q_k: DMatrix<f64>, estimate: KfEstimate) {
        let (f_rows, f_cols) = (f_k.nrows(), f_k.ncols());
        let (q_rows, q_cols) = (q_k.nrows(), q_k.ncols());

        assert_eq!(f_rows, f_cols, "dimensions: F is not square");
        assert_eq!(f_rows, q_cols, "invalid F/Q dimensions");
        assert_eq!(f_cols, q_rows, "invalid F/Q dimensions");

        // prediction
        let x_k = f_k.clone() * estimate.x;
        let f_k_t = f_k.transpose();
        let p_k = f_k * estimate.p * f_k_t + q_k;

        self.predicted = KfEstimate::new(&x_k, &p_k);
        self.initialized = true;
    }

    /// Reset this [Kalman] filter
    pub fn reset(&mut self) {
        self.initialized = false;
        self.predicted.reset();
    }

    /// Initializes this filter as an error-state filter, from the
    /// covariance of the initial (least squares) estimate.
    ///
    /// The error-state formulation keeps the nominal state outside of the filter:
    /// after each update, the correction is applied to the nominal state
    /// and the error-state is reset to zero. Only the covariance lives here.
    pub fn initialize_error_state(&mut self, p_0: &DMatrix<f64>) {
        let n = p_0.nrows();
        self.predicted = KfEstimate {
            x: DVector::<f64>::zeros(n),
            p: p_0.clone(),
        };
        self.initialized = true;
    }

    /// Current covariance
    pub fn covariance(&self) -> &DMatrix<f64> {
        &self.predicted.p
    }

    /// Error-state prediction: P⁻ = F P Fᵀ + Q (the predicted error is null).
    pub fn predict_error_state(&mut self, f_k: &DMatrix<f64>, q_k: &DMatrix<f64>) {
        assert!(self.initialized, "internal error: filter not initialized!");
        let p = f_k * &self.predicted.p * f_k.transpose() + q_k;
        let n = p.nrows();
        self.predicted = KfEstimate {
            x: DVector::<f64>::zeros(n),
            p,
        };
    }

    /// Error-state measurement update (standard form + Joseph covariance update).
    ///
    /// ## Input
    /// - g_k: observation (design) matrix H
    /// - w_k: weight matrix, inverse of the measurement covariance R (diagonal)
    /// - z_k: innovation vector (observed - computed, at the nominal state)
    ///
    /// ## Output
    /// - [KfEstimate] where x is the state correction (δx) and P the posterior covariance.
    ///   The internal error-state is reset (δx applied by the caller).
    pub fn update_error_state(
        &mut self,
        g_k: &DMatrix<f64>,
        w_k: &DMatrix<f64>,
        z_k: &DVector<f64>,
    ) -> Result<KfEstimate, Error> {
        assert!(self.initialized, "internal error: filter not initialized!");

        let n = self.predicted.p.nrows();
        let m = z_k.nrows();

        if g_k.ncols() != n || g_k.nrows() != m || w_k.nrows() != m || w_k.ncols() != m {
            return Err(Error::MatrixDimension);
        }

        // R = W⁻¹ (diagonal weights)
        let mut r_k = DMatrix::<f64>::zeros(m, m);

        for i in 0..m {
            let w = w_k[(i, i)];
            if w <= 0.0 {
                return Err(Error::MatrixFormationError);
            }
            r_k[(i, i)] = 1.0 / w;
        }

        let p_minus = &self.predicted.p;
        let gt = g_k.transpose();

        let s_k = g_k * p_minus * &gt + &r_k;
        let s_inv = s_k.try_inverse().ok_or(Error::MatrixInversion)?;

        let k_k = p_minus * &gt * s_inv;

        let dx = &k_k * z_k;

        // Joseph form: numerically stable, keeps P symmetric positive
        let i_kh = DMatrix::<f64>::identity(n, n) - &k_k * g_k;
        let mut p_plus = &i_kh * p_minus * i_kh.transpose() + &k_k * r_k * k_k.transpose();

        // enforce symmetry
        p_plus = (&p_plus + p_plus.transpose()) * 0.5;

        self.predicted = KfEstimate {
            x: DVector::<f64>::zeros(n),
            p: p_plus.clone(),
        };

        Ok(KfEstimate { x: dx, p: p_plus })
    }

    /// Run this [Kalman] filter, returning new [KfEstimate].
    ///
    /// ## Input
    /// - f_k: Dynamics [DMatrix]
    /// - g_k: G [DMatrix]
    /// - w_k: W [DMatrix]
    /// - q_k: Q [DMatrix]
    /// - y_k: Measurement [DVector]
    pub fn run(
        &mut self,
        f_k: &DMatrix<f64>,
        g_k: &DMatrix<f64>,
        w_k: &DMatrix<f64>,
        q_k: &DMatrix<f64>,
        y_k: &DVector<f64>,
    ) -> Result<KfEstimate, Error> {
        let (w_rows, w_cols) = (w_k.nrows(), w_k.ncols());
        let (g_rows, g_cols) = (g_k.nrows(), g_k.ncols());
        let (f_rows, f_cols) = (f_k.nrows(), f_k.ncols());
        let (q_rows, q_cols) = (q_k.nrows(), q_k.ncols());

        let y_rows = y_k.nrows();

        if !self.initialized {
            panic!("internal error: filter not initialized!");
        }

        assert_eq!(w_rows, w_cols, "W is not square");
        assert_eq!(f_rows, f_cols, "F is not square");

        assert_eq!(f_rows, q_cols, "invalid F/Q dimensions");
        assert_eq!(f_cols, q_rows, "invalid F/Q dimensions");

        assert_eq!(q_cols, g_cols, "invalid Q/G dimensions!");
        assert_eq!(f_cols, g_cols, "invalid F/G dimensions!");
        assert_eq!(y_rows, g_rows, "invalid Y/G dimensions!");
        assert_eq!(y_rows, w_rows, "invalid Y/W dimensions!");

        let gt = g_k.transpose();

        let p_inv = self
            .predicted
            .p
            .clone()
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        let p_k = gt.clone() * w_k.clone();
        let p_k = p_k * g_k;
        let p_k = p_k + p_inv.clone();
        let p_k = p_k.try_inverse().ok_or(Error::MatrixInversion)?;

        let p_inv_x = p_inv.clone() * self.predicted.x.clone();

        let x_k = gt * w_k;
        let x_k = x_k * y_k;
        let x_k = x_k + p_inv_x;
        let x_k = p_k.clone() * x_k;

        // prediction
        let x_k1 = f_k.clone() * x_k.clone();

        let p_k1 = f_k.clone() * p_k.clone() * f_k.transpose() + q_k;

        let x_k1 = DVector::<f64>::from_row_slice(x_k1.as_slice());
        let p_k1 = DMatrix::<f64>::from_column_slice(x_k1.nrows(), p_k1.ncols(), p_k1.as_slice());

        self.predicted = KfEstimate { x: x_k1, p: p_k1 };

        Ok(KfEstimate { x: x_k, p: p_k })
    }

    /// Run this [Kalman] filter, returning new [KfEstimate].
    ///
    /// ## Input
    /// - f_k: Dynamics [DMatrix]
    /// - g_k: G [DMatrix]
    /// - w_k: W [DMatrix]
    /// - q_k: Q [DMatrix]
    /// - y_k: Measurement [OVector]
    pub fn run_static<D: DimName>(
        &mut self,
        f_k: &OMatrix<f64, D, D>,
        g_k: &OMatrix<f64, D, D>,
        w_k: &OMatrix<f64, D, D>,
        q_k: &OMatrix<f64, D, D>,
        y_k: &OVector<f64, D>,
    ) -> Result<KfEstimate, Error>
    where
        DefaultAllocator: nalgebra::allocator::Allocator<D>,
        DefaultAllocator: nalgebra::allocator::Allocator<D, D>,
    {
        if !self.initialized {
            panic!("internal error: filter not initialized!");
        }

        let gt = g_k.transpose();

        let p_inv = self
            .predicted
            .p
            .clone()
            .try_inverse()
            .ok_or(Error::MatrixInversion)?;

        let p_k = gt.clone() * w_k.clone();
        let p_k = p_k * g_k;
        let p_k = p_k + p_inv.clone();
        let p_k = p_k.try_inverse().ok_or(Error::MatrixInversion)?;

        let p_inv_x = p_inv.clone() * self.predicted.x.clone();

        let x_k = gt * w_k;
        let x_k = x_k * y_k;
        let x_k = x_k + p_inv_x;
        let x_k = p_k.clone() * x_k;

        // prediction
        let x_k1 = f_k.clone() * x_k.clone();
        let p_k1 = f_k.clone() * p_k.clone() * f_k.transpose() + q_k;

        self.predicted = KfEstimate::from_static(x_k1, p_k1);

        Ok(KfEstimate::from_static(x_k, p_k))
    }
}

#[cfg(test)]
mod test {
    use super::KfEstimate;

    use nalgebra::{DimName, OMatrix, OVector, U6};

    #[test]
    fn kf_estimate_from_static() {
        let x = OVector::<f64, U6>::from_column_slice(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        let p = OMatrix::<f64, U6, U6>::from_column_slice(&[
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, // 1
            7.0, 8.0, 9.0, 10.0, 11.0, 12.0, // 2
            13.0, 14.0, 15.0, 16.0, 17.0, 18.0, // 3
            19.0, 20.0, 21.0, 22.0, 23.0, 24.0, // 4
            25.0, 26.0, 27.0, 28.0, 29.0, 30.0, // 5
            31.0, 32.0, 33.0, 34.0, 35.0, 36.0, // 6
        ]);

        let kfe = KfEstimate::from_static(x, p);

        for i in 0..U6::DIM {
            let x = kfe.x[i] as usize;
            assert_eq!(x, i + 1);

            for j in 0..U6::DIM {
                let p = kfe.p[(i, j)] as usize;
                assert_eq!(p, i + 1 + j * U6::DIM);
            }
        }
    }
}
